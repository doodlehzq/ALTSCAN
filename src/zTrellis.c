/* -*- Mode: C; tab-width: 4; c-basic-offset: 4; indent-tabs-mode: t -*- */  
/******************************************************************************\
  zTrellis.c - part of the ZOE library for genomic analysis
  
  Copyright (C) 2001-2002 Ian F. Korf
  Modified by Zhiqiang Hu in 2012
\******************************************************************************/

#ifndef ZOE_TRELLIS_C
#define ZOE_TRELLIS_C

#include <assert.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "zTrellis.h"
#include "zFastaFile.h"

/* static const int PADDING = 50; */
static const score_t MIN_INIT_SCORE = -10000;
extern coor_t any_new_cpoint;


void zShowTrellis(zTrellis *trellis);
score_t zScoreInternalStateRange(zTrellis *trellis, int state, coor_t start, coor_t end){
	zHMM* hmm = trellis->hmm;
	score_t score = 0;
	score_t tscore;
	coor_t i;

	if(hmm->state[state].type != INTERNAL &&
	   hmm->state[state].type != GINTERNAL){
		zDie("State %s is not an internal state\n",zStrIdx2Char(trellis->hmm->state[state].name));
	}
	
	if(hmm->state[state].type == INTERNAL){
		tscore = zGetTransitionScore(trellis->hmm,state,state,trellis->tiso_group);
	}
	else{
		tscore = zGetIntergenicContinueScore(trellis);
	}
	
	score = zScoreInternalState(trellis,state,start,true);
	
	for(i = start+1; i <= end; i++) {
		score += tscore;
		score += zScoreInternalState(trellis,state,i,false);
	}
	
	return score;
}


void zTracePartialTrellis (zTrellis *trellis, coor_t start, coor_t end, int state, zSFVec **sfv,score_t *max_score){
	zHMM      *hmm = trellis->hmm;
	int        max_state, trace;
	int        top_trace_top[top_num];
	coor_t     i;
	zSfeature  istate;
	score_t    max_score_temp[top_num];
	int        m,n1,n2;
	int        top_state[top_num];
	int        top_state_temp[top_num];
	int        top_trace_top_temp[top_num];
	int        trace_top;
	bool       flag_tr;
	if (state < 0) {
		if(!zOption("top")){
			max_score[0] = MIN_SCORE;
			for (max_state = 0; max_state < hmm->states; max_state++) {
				if (trellis->cell[max_state] == NULL) continue;
		
				max_score_temp[0] = trellis->cell[max_state][end].score[0];
				if (max_score_temp[0] > max_score[0]) {
					state = max_state; max_score[0]= max_score_temp[0];
				}
			}
		}
		else{
			for(m=0;m<top_num;m++){
				max_score[m] = MIN_SCORE;
				top_state[m] = -1;
				top_trace_top[m] = -1;
			}
			for (max_state = 0; max_state < hmm->states; max_state++) {
				if (trellis->cell[max_state] == NULL) continue;
				for(m=0;m<top_num;m++){
					max_score_temp[m] = max_score[m]; 
					top_state_temp[m] = top_state[m];
					top_trace_top_temp[m] = top_trace_top[m];
				}
				n1=0;n2=0;
				for(m=0;m<top_num;m++){
					if(trellis->cell[max_state][end].score[n1] != MIN_SCORE){
						flag_tr = (max_score_temp[n2] == MIN_SCORE);
						if(!flag_tr)
							flag_tr = (max_score_temp[n2] < 
									   trellis->cell[max_state][end].score[n1]);
						if (flag_tr) {
							top_state[m] = max_state; 
							max_score[m] = trellis->cell[max_state][end].score[n1];
							top_trace_top[m] = n1;
							n1++;
						}
						else{
							top_state[m] = top_state_temp[n2];
							max_score[m] = max_score_temp[n2];
							top_trace_top[m] = top_trace_top_temp[n2];
							n2++;
						}
					}
					else{
					   	top_state[m] = top_state_temp[n2];
					   	max_score[m] = max_score_temp[n2];
					   	top_trace_top[m] = top_trace_top_temp[n2];
					   	n2++;
					}
				}
			}
		}
	}
	/*zShowTrellis(trellis);
	zDie("show trellis now!\n");*/
	if(!zOption("top")){
		i = end;
		if (-1 == trellis->cell[state][i].trace[0]) {
			return; /* the sequence can't end in this state */
		}
	
		while (i > start) {
			trace = trellis->cell[state][i].trace[0];
			if (-1 == trace) {
				zDie("Impossible to perform traceback from this state.");
			}

			istate.name   = hmm->state[state].name;
			istate.start  = i - trellis->cell[state][i].length[0] + 1;
			istate.end    = i;
			istate.strand = hmm->state[state].strand;

			zChar2FragFrame(trellis->cell[state][i].frag_frame[0], &istate);

			istate.state  = state;
			istate.group  = NULL;

			istate.score = trellis->cell[state][istate.end].score[0]
				- trellis->cell[trace][istate.start-1].score[0]
				- zGetTransitionScore(trellis->hmm,trace,state,trellis->tiso_group);
			/*EVAN now removing transition scores from output */

			zPushSFVec(sfv[0], &istate);

			assert(trellis->cell[state][i].length[0] >  0);
			i -= trellis->cell[state][i].length[0];
			state=trace;
		}
	}
	else{
		int ttt;

		for(m=0;m<top_num;m++){
			if(max_score[m] != MIN_SCORE){
				if(top_trace_top[m]<0 || top_trace_top[m] >= top_num){
					zDie("top_trace_top out of range!\n");
				}
			}
		}
		
		for(m=0;m<top_num;m++){
			if(max_score[m] == MIN_SCORE) break;
			i = end;
			
			/*	if (-1 == trellis->cell[state][i].trace[trace_top]) {
			return; 
			}*/ 
			state = top_state[m];
			ttt = top_trace_top[m];

			while (i > start) {
			
			   	trace = trellis->cell[state][i].trace[ttt];
				trace_top = trellis->cell[state][i].trace_top[ttt];
				if(trace_top == -1) zDie("trace_top error!\n");

				istate.name   = hmm->state[state].name;
				istate.start  = i - trellis->cell[state][i].length[ttt] + 1;
				istate.end    = i;
				istate.strand = hmm->state[state].strand;

				zChar2FragFrame(trellis->cell[state][i].frag_frame[ttt], &istate);

				istate.state  = state;
				istate.group  = NULL;

			   	istate.score = trellis->cell[state][istate.end].score[ttt]
					- trellis->cell[trace][istate.start-1].score[trace_top]
					- zGetTransitionScore(trellis->hmm,trace,state,trellis->tiso_group);
				/*EVAN now removing transition scores from output */
				zPushSFVec(sfv[m], &istate);

				assert(trellis->cell[state][i].length[ttt] >  0);
				i -= trellis->cell[state][i].length[ttt];
				state = trace;
			   	ttt = trace_top;
			}
		}
	}
}

void zTraceTrellis (zTrellis *trellis, int state, zSFVec **sfv, score_t *score) {
    zTracePartialTrellis(trellis, PADDING, trellis->dna->length - 1 - PADDING, state, sfv, score);
}

static void zAllocViterbiVars(zTrellis *trellis) {
	coor_t    j;
	int       i;
	int       k;
	size_t size;
	if (NULL != trellis->cell) return;
	trellis->cell = zCalloc(trellis->hmm->states, sizeof(zTrellisCell*), 
							"zAllocViterbi cells");

	/* Allocating backward link arrays for all states                 */
	/* But only those referring to external states will be used later */
	trellis->extpos = zMalloc((trellis->hmm->states*sizeof(zIVec)), "zAllocViterbi extpos");

	size = trellis->dna->length * sizeof(zTrellisCell);
	for (i = 0; i < trellis->hmm->states; i++) {
		trellis->cell[i] = zMalloc(size, "zAllocViterbi cells[i]");
		zInitIVec(&trellis->extpos[i], 1);

		for (j = 0; j < trellis->dna->length; j++){
			for(k=0; k < top_num; k++){
				trellis->cell[i][j].score[k] = MIN_SCORE;
			}
		}
	}
}
	
void zFreeViterbiVars(zTrellis *trellis) {
	int    i;

	if (NULL == trellis->cell) return;
	for (i = 0; i < trellis->hmm->states; i++) {
		if (NULL == trellis->cell[i]) continue;
		zFree(trellis->cell[i]);
		trellis->cell[i] = NULL;

		zFreeIVec(&trellis->extpos[i]);
	}

	zFree(trellis->cell);
	trellis->cell = NULL;

	zFree(trellis->extpos);
	trellis->extpos = NULL;
}

void zAllocForwardVars(zTrellis *trellis) {
	int i;
	if (NULL != trellis->forward) return;
	trellis->forward = zCalloc(trellis->hmm->states, sizeof(score_t*), 
							   "zAllocForward forward");
	for (i = 0; i < trellis->hmm->states; i++) {
		trellis->forward[i] = zMalloc(trellis->dna->length * sizeof(score_t),
									  "zAllocForward forward[i]");
	}
}

void zFreeForwardVars(zTrellis *trellis) {
	int i;
	if (NULL == trellis->forward) return;
	for (i = 0; i < trellis->hmm->states; i++) {
		if (NULL == trellis->forward[i]) continue;
		zFree(trellis->forward[i]);
		trellis->forward[i] = NULL;
	}
	zFree(trellis->forward);
	trellis->forward = NULL;
}

static void zAllocBackwardVars(zTrellis *trellis) {
	int i;
	if (NULL != trellis->forward) return;
	trellis->backward = zCalloc(trellis->hmm->states, sizeof(score_t*),
								"zAllocBackwardVars backward");
	for (i = 0; i < trellis->hmm->states; i++) {
		trellis->backward[i] = zMalloc(trellis->dna->length * sizeof(score_t),
									   "zAllocBackward forward[i]");
	}
}

void zFreeBackwardVars(zTrellis *trellis) {
	int i;
	if (NULL == trellis->backward) return;
	for (i = 0; i < trellis->hmm->states; i++) {
		if (NULL == trellis->backward[i]) continue;
		zFree(trellis->backward[i]);
		trellis->backward[i] = NULL;
	}
	zFree(trellis->forward);
	trellis->forward = NULL;
}

static void zAllocFactories(zTrellis *trellis, bool conseq_enabled) {
	zFeatureFactory **ffactory;
	zHMM_State        *state;
	zIVec           **fmap5,**fmap3,*jumps;

	zHMM    *hmm = trellis->hmm;
	int      i,j,k,insert;
	int      ns;

	/* create one zStopSeq for each strand */
	trellis->stopseq = zMalloc(sizeof(zStopSeq),"zAllocFactories stopseq");
	trellis->rstopseq = zMalloc(sizeof(zStopSeq),"zAllocFactories rstopseq");
	zInitStopSeq(trellis->stopseq,trellis->unmasked_dna == NULL ? trellis->dna :
				 trellis->unmasked_dna);
	zInitStopSeq(trellis->rstopseq,trellis->unmasked_rdna == NULL ? trellis->rdna :
				 trellis->unmasked_rdna);
	
	
	/* create factories for non-internal states */
	trellis->factory = zCalloc(hmm->feature_count, sizeof(zFeatureFactory*),
							   "zAllocFactories: factory[]");
	trellis->rfactory = zCalloc(hmm->feature_count, sizeof(zFeatureFactory*),
								"zAllocFactories: ractory[]");
	trellis->fmap5   = zCalloc(hmm->feature_count, sizeof(zIVec*),
														 "zAllocFactories: fmap5[]");
	trellis->rfmap5  = zCalloc(hmm->feature_count, sizeof(zIVec*),
														 "zAllocFactories: rfmap5[]");
	trellis->fmap3   = zCalloc(hmm->feature_count, sizeof(zIVec*),
														 "zAllocFactories: fmap3[]");
	trellis->rfmap3  = zCalloc(hmm->feature_count, sizeof(zIVec*),
														 "zAllocFactories: rfmap3[]");
	for (i = 0; i <  hmm->states; i++) {
		state = &hmm->state[i];
		if (NULL == state->ffactory) continue;

		ns = ('-' == state->strand);
		ffactory = (ns) ? &trellis->rfactory[state->model] : 
			&trellis->factory[state->model];
		fmap5 = (ns) ? &trellis->rfmap5[state->model] : &trellis->fmap5[state->model];
		fmap3 = (ns) ? &trellis->rfmap3[state->model] : &trellis->fmap3[state->model];
		if (NULL != *ffactory){
			zPushIVec(*fmap3,i);
			jumps = trellis->hmm->jmap[i];
			for(j = 0;j < jumps->size; j++){
				insert = 1;
				for(k = 0;k < (*fmap5)->size; k++){
					if((*fmap5)->elem[k] == jumps->elem[j]){
						insert = 0;
						break;
					}
				}
				if(insert == 1){
					zPushIVec(*fmap5,jumps->elem[j]);
				}
			}
			continue;
		}
		*fmap5 = zMalloc(sizeof(zIVec),"zAllocFacotries: fmap5");
		*fmap3 = zMalloc(sizeof(zIVec),"zAllocFacotries: fmap3");
		zInitIVec(*fmap5,3);
		zInitIVec(*fmap3,3);
		zPushIVec(*fmap3,i);
		jumps = trellis->hmm->jmap[i];
		for(j = 0;j < jumps->size; j++){
			zPushIVec(*fmap5,jumps->elem[j]);
		}
		*ffactory = zMalloc(sizeof(zFeatureFactory),"zAllocFac: factory[i]");
		state->ffactory(*ffactory,
						(ns) ? trellis->rdna     : trellis->dna, 
						(ns) ? trellis->rscanner : trellis->scanner, 
						(ns) ? trellis->rconseqscanner: trellis->conseqscanner,
							hmm->feature_count, 
						state->strand, 
						conseq_enabled,
							(ns) ? trellis->rstopseq : trellis->stopseq
						);
	
	}
	/* allcate external arrays */
	trellis->external = zMalloc(hmm->feature_count*sizeof(zSFList*),
								"zAllocFactories: external[]");
	trellis->rexternal = zMalloc(hmm->feature_count*sizeof(zSFList*),
								 "zAllocFactories: rexternal[]");
	for (i = 0; i <  hmm->feature_count; i++) {
		if(trellis->factory[i] != NULL){
			trellis->external[i] = zMalloc(sizeof(zSFList),"zAllocFactories: external[i]");
			zInitSFList(trellis->external[i]);
		}
		if(trellis->rfactory[i] != NULL){
			trellis->rexternal[i] = zMalloc(sizeof(zSFList),"zAllocFactories: rexternal[i]");
			zInitSFList(trellis->rexternal[i]);
		}
	}
}

static void zFreeFactories(zTrellis *trellis) {
	int i;

	for (i = 0; i < trellis->hmm->feature_count; i++) {
		if (NULL != trellis->factory[i]) {
			zFreeFeatureFactory(trellis->factory[i]);
			zFree(trellis->factory[i]);
			zFreeSFList(trellis->external[i]);
			zFree(trellis->external[i]);
			zFreeIVec(trellis->fmap5[i]);
			zFreeIVec(trellis->fmap3[i]);
			zFree(trellis->fmap5[i]);
			zFree(trellis->fmap3[i]);
		}
		if (NULL != trellis->rfactory[i]) {
			zFreeFeatureFactory(trellis->rfactory[i]);
			zFree(trellis->rfactory[i]);
			zFreeSFList(trellis->rexternal[i]);
			zFree(trellis->rexternal[i]);
			zFreeIVec(trellis->rfmap5[i]);
			zFreeIVec(trellis->rfmap3[i]);
			zFree(trellis->rfmap5[i]);
			zFree(trellis->rfmap3[i]);
		}
	}
	zFree(trellis->fmap5);
	zFree(trellis->rfmap5);
	zFree(trellis->fmap3);
	zFree(trellis->rfmap3);
	zFree(trellis->factory);
	zFree(trellis->rfactory);
	zFree(trellis->external);
	zFree(trellis->rexternal);
	trellis->external = trellis->rexternal = NULL;
	trellis->factory = trellis->rfactory = NULL;

	zFreeStopSeq(trellis->stopseq);
	zFreeStopSeq(trellis->rstopseq);
	zFree(trellis->stopseq);
	zFree(trellis->rstopseq);
	trellis->stopseq = NULL;
	trellis->rstopseq = NULL;
}

void zInitTrellis (zTrellis *trellis, char* dna_file, zHMM *hmm, char* conseq_file, 
				   char* unmasked_dna_file, char* snp_file, char* p_top_num) {
	int         i;
	zDNA       *dna;
	bool        conseq_enabled;

	/* clear out pointers */
	trellis->dna       = NULL;
	trellis->rdna      = NULL;
	trellis->conseq    = NULL;
	trellis->rconseq    = NULL;
	trellis->hmm       = NULL;
	trellis->fexternal = NULL;
	trellis->cell      = NULL;
	trellis->forward   = NULL;
	trellis->backward  = NULL;
	trellis->scanner   = NULL;
	trellis->factory   = NULL;
	trellis->extpos    = NULL;
	trellis->external  = NULL;
	trellis->rexternal = NULL;

	if(conseq_file != NULL) {
		conseq_enabled = true;
	}
	else {
		conseq_enabled = false;
	}
	
	/* initial setup */
	trellis->dna  = zMalloc(sizeof(zDNA), "zInitTrellis dna");
	trellis->rdna = zMalloc(sizeof(zDNA), "zInitTrellis rdna");
	zInitDNA(trellis->dna);
	zInitDNA(trellis->rdna);
	dna = trellis->dna;
	zLoadDNAFromFasta(dna,dna_file,snp_file);
	zSetDNAPadding(dna,PADDING);
	trellis->padding = PADDING;
	zCopyDNA(trellis->dna, trellis->rdna);
	zAntiDNA(trellis->rdna);
	if(unmasked_dna_file != NULL){
		trellis->unmasked_dna  = zMalloc(sizeof(zDNA), "zInitTrellis unmasked_dna");
		trellis->unmasked_rdna = zMalloc(sizeof(zDNA), "zInitTrellis unmasked_rdna");
		zInitDNA(trellis->unmasked_dna);
		zInitDNA(trellis->unmasked_rdna);
		dna = trellis->dna;
		zLoadDNAFromFasta(trellis->unmasked_dna,unmasked_dna_file,snp_file);
		zSetDNAPadding(trellis->unmasked_dna,PADDING);
		zCopyDNA(trellis->unmasked_dna, trellis->unmasked_rdna);
		zAntiDNA(trellis->unmasked_rdna);
	}
	else{
		trellis->unmasked_dna = NULL;
		trellis->unmasked_rdna = NULL;
	}
	if(p_top_num!=NULL){
	}
	if(conseq_enabled) {
		trellis->conseq = zMalloc(sizeof(zConseq), "zInitTrellis conseq");
		trellis->rconseq = zMalloc(sizeof(zConseq), "zInitTrellis rconseq");
		zInitConseq(trellis->conseq,3);
		zInitConseq(trellis->rconseq,3);
		zLoadConseqFromFasta(trellis->conseq,conseq_file);
		zSetConseqPadding(trellis->conseq,PADDING);
		zCopyConseq(trellis->conseq, trellis->rconseq);
		zReverseConseq(trellis->rconseq);
	}
	
	trellis->hmm  = hmm;

	/* isocore group */
	trellis->tiso_group = -1;
	if(trellis->hmm->iso_transitions > 0) {
		for(i=0;i<trellis->hmm->iso_transitions;i++){
			if(zGetDNAGC(trellis->dna) < trellis->hmm->iso_transition[i]){
				trellis->tiso_group = i;
				break;
			}
		}
	}
	trellis->iiso_group = -1;
	if(hmm->iso_states > 0){
		for(i = 0; i < hmm->iso_states; i++){
			if(zGetDNAGC(trellis->dna) < hmm->iso_state[i]){
				trellis->iiso_group=i;
				break;
			}
		}
	}	

	/* create scanners */
	trellis->scanner  = zCalloc(hmm->feature_count, sizeof(zScanner*), 
								"zInitTrellis: scanner[]");
	trellis->rscanner = zCalloc(hmm->feature_count, sizeof(zScanner*), 
								"zInitTrellis: rscanner[]");
	
	for (i = 0; i < hmm->feature_count; i++) {
		if (hmm->mmap[i] == NULL) continue;
		if (hmm->mmap[i]->seq_type == DNA) {
			trellis->scanner[i] = zMalloc(sizeof(zScanner), "zInitTrellis scan");
			zInitScanner(trellis->scanner[i], trellis->dna->seq, hmm->mmap[i]);
			trellis->rscanner[i] = zMalloc(sizeof(zScanner), "zInitTrellis rscan");
			zInitScanner(trellis->rscanner[i], trellis->rdna->seq, hmm->mmap[i]);

			/* Pre-compute DNA scanners for EXPLICIT    */
			/* states as a "running sum". Required by   */
			/* the implementation of EXPLICIT states    */
			/* (see zScoreScanners in zTransition.c for */
			/* details).                                */
			
			if (zUsedInExplicitState(hmm, i)){
				zPreComputeScanner(trellis->scanner[i]);
				zPreComputeScanner(trellis->rscanner[i]);
			}
		}
		else {
			zDie("non-DNA model in sequence models\n");
		}
	}

	/* create conseq scanners */
	if(conseq_enabled) {
		trellis->conseqscanner = 
			zCalloc(hmm->feature_count, sizeof(zScanner*), 
					"zInitTrellis: conseqscanner[]");
		trellis->rconseqscanner = 
			zCalloc(hmm->feature_count, sizeof(zScanner*), 
					"zInitTrellis: rconseqscanner[]");
	} else {
		trellis->conseqscanner = trellis->rconseqscanner = NULL;
	}
	for (i = 0; i < hmm->feature_count; i++) {
		if (NULL == hmm->cmmap[i]) continue;
		if (hmm->cmmap[i]->seq_type == CONSEQ) {
			if (conseq_enabled) {
				trellis->conseqscanner[i] = zMalloc(sizeof(zScanner), 
													"zInitTrellis conseqscan");
				zInitScanner(trellis->conseqscanner[i], trellis->conseq->seq, 
							 hmm->cmmap[i]);
				trellis->rconseqscanner[i] = zMalloc(sizeof(zScanner), 
													 "zInitTrellis rcscan");
				zInitScanner(trellis->rconseqscanner[i], 
							 trellis->rconseq->seq, hmm->cmmap[i]);

                                /* Pre-compute conseq scanners for EXPLICIT */
                                /* states as a "running sum". Required by   */
                                /* the implementation of EXPLICIT states    */
                                /* (see zScoreScanners in zTransition.c for */
                                /* details).                                */

				/*if (zUsedInExplicitState(hmm, i)){ EVAN removed this for speed*/
				zPreComputeScanner(trellis->conseqscanner[i]);
				zPreComputeScanner(trellis->rconseqscanner[i]);
				/*}*/

			}
		} else {
			zDie("non-conservation model in list of conservation models");
		}
	}

	zAllocFactories(trellis, conseq_enabled);
}

void zFreeTrellis (zTrellis *trellis) {
	int i;

	
	for (i = 0; i < trellis->hmm->feature_count; i++) {
		if (trellis->scanner[i]  != NULL) {
			zFreeScanner(trellis->scanner[i]);
			zFree(trellis->scanner[i]);
		}
		if (trellis->rscanner[i] != NULL) {
			zFreeScanner(trellis->rscanner[i]);
			zFree(trellis->rscanner[i]);
		}

		if (trellis->conseq) {
			if (trellis->conseqscanner[i] != NULL) {
				zFreeScanner(trellis->conseqscanner[i]);
				zFree(trellis->conseqscanner[i]);
			}
			if (trellis->rconseqscanner[i] != NULL) {
				zFreeScanner(trellis->rconseqscanner[i]);
				zFree(trellis->rconseqscanner[i]);
			}
		}
	}

	zFree(trellis->scanner);        trellis->scanner = NULL; 
	zFree(trellis->rscanner);       trellis->rscanner = NULL;
	zFree(trellis->conseqscanner);  trellis->conseqscanner = NULL;
	zFree(trellis->rconseqscanner); trellis->rconseqscanner = NULL;

	
	zFreeFactories(trellis);

	zFreeViterbiVars(trellis);
	zFreeForwardVars(trellis);
	zFreeBackwardVars(trellis);
	
	zFreeDNA(trellis->rdna);
	zFree(trellis->rdna);
	trellis->rdna = NULL;
	zFreeDNA(trellis->dna);
	zFree(trellis->dna);
	trellis->dna = NULL;

	if (trellis->conseq != NULL) {
		zFreeConseq(trellis->conseq);
		zFreeConseq(trellis->rconseq);
		zFree(trellis->conseq);
		zFree(trellis->rconseq);
	}
}

void zRunPartialViterbiAndForward(zTrellis *trellis, coor_t start, coor_t end){
	zHMM         *hmm = trellis->hmm;

	coor_t        i;        /* iterator for sequence */
	int           j;        /* iterator for internal states */
	int           k;        /* iterator for previous states */
	zIVec*        jumps;

	/* induction */
	zTrace2("beginning induction\n");
	for (i = start; i <= end; i++) {
		for (j = 0; j < hmm->states; j++) {
			if (NULL == trellis->cell[j]) continue;
			/* 			if (NULL == trellis->forward[j]) continue; */

			trellis->cell[j][i].score[0] = MIN_SCORE;
			trellis->cell[j][i].trace[0] = -1;
			trellis->cell[j][i].length[0] = 0;
			/* 			trellis->forward[j][i] = MIN_SCORE; */

			jumps = hmm->jmap[j];
			for(k = 0; k < jumps->size; k++) {
				zGetTransFunc(hmm->state[j].type)
					(trellis, jumps->elem[k], j, i);
			}
		}
	}
}

zSFVec** zRunViterbiAndForward (zTrellis *trellis, score_t* path_score) {
	zHMM         *hmm = trellis->hmm;
	zDNA         *dna = trellis->dna;

	coor_t        i,l,tmp_start;
	int           j;        /* iterator for internal states */
	int           k;        /* iterator for previous states */
	int           m;
	/*	int           m;    */    /*top mode*/
	int           g;
	int           *nrange, iso_group, u;
	zDurationGroup *group;
	zIVec*        jumps;
	zSFVec        **sfv;
	zSFList       *sfl;
	zSfeature     *f;
	size_t        size;
	/*-----------------------------------------------
	  |            [0] [1] [2] [3] [4] [5] [6] [7] ...
	  |(None)  [0]
	  |(Misc)  [1]              X
	  |(Inter) [2]
	  |:                  X = trellis->cell[1][3] 
	  |:                               cell[j][i]
	  ------------------------------------------------*/
	 	 
	/* 	Viterbi and Forward Alg initialization */
	zTrace2("initializing Viterbi vars\n");

	zAllocViterbiVars(trellis);

	sfl = zMalloc(sizeof(zSFList),"zRunViterbiAndForward: sfl");
	zInitSFList(sfl);

	/* EVAN what is this nrange crap?*/
	nrange = zMalloc((hmm->states*sizeof(int)), "ZRunViterbiAndForward nrange");
	for (j = 0; j < hmm->states; j++){
		group = hmm->dmap[hmm->state[j].duration];
		iso_group = zGetDurationIsochoreGroup(group, zGetDNAGC(dna));
		
		if ((hmm->state[j].type == EXPLICIT) 
			|| (hmm->state[j].type == INTERNAL)){
			nrange[j] = group->duration[iso_group].distributions;
		}
	}

	/* EVAN there must be a better way to do this explicit stuff */

 	zTrace2("starting initilization\n");
	for (j = 0; j < hmm->states; j++) {
		for (i = 0; i < PADDING; i++) {
			for(m=0;m<top_num;m++){
				if(m==0){
					trellis->cell[j][i].score[m] = zGetInitProb(trellis->hmm, j, trellis->iiso_group);
					trellis->cell[j][i].trace[m] = j;
					trellis->cell[j][i].length[m] = 0;
					trellis->cell[j][i].trace_top[m] = 0;
				}
				else{
					trellis->cell[j][i].score[m] = MIN_SCORE;
					trellis->cell[j][i].trace[m] = -1;
					trellis->cell[j][i].length[m] = -1;
					trellis->cell[j][i].trace_top[m] = -1;
				}
			}
			/* Allows introns to extend to the start of the */
			/* sequence when modeled as an EXPLICIT state   */
			if ((trellis->cell[j][i].score[0] > MIN_SCORE)
				&& (hmm->state[j].type == EXTERNAL)
				&& ((i+1) == PADDING)){
				zPushIVec(&trellis->extpos[j], i);
			}
			
			if (trellis->cell[j][i].score[0] == MIN_SCORE) {
				trellis->cell[j][i].score[0] = MIN_INIT_SCORE;
			}
			
			if (((i+1) == PADDING)
				&& ((hmm->state[j].type == EXPLICIT)
					|| (hmm->state[j].type == INTERNAL))){
				trellis->cell[j][i].submax = zMalloc(sizeof(zSFVec),
													 "zRunViterbiAndForward submax");
				zInitSFVec(trellis->cell[j][i].submax, nrange[j]);
				trellis->cell[j][i].submax->size = nrange[j];
				trellis->cell[j][i].submax->last =
					&trellis->cell[j][i].submax->elem[nrange[j]-1];
				
				for (u = 0; u < nrange[j]; u++){
					trellis->cell[j][i].submax->elem[u].score =
						trellis->cell[j][i].score[0];
					trellis->cell[j][i].submax->elem[u].intrinsic_score =
						trellis->cell[j][i].score[0];
					trellis->cell[j][i].submax->elem[u].from_state = -1;
					trellis->cell[j][i].submax->elem[u].end = i;
					trellis->cell[j][i].submax->elem[u].start = i;
				}
			}
		}
	}
	
	/* induction - first pass */
	zTrace2("beginning first induction\n");
	for (i = PADDING; i < dna->length; i++) {
		/* create external links ending at i */
		for (j = 0; j < trellis->hmm->feature_count; j++) {
			l = i;
			if(trellis->factory[j] != NULL){
				zResetSFList(trellis->external[j]);		
				trellis->factory[j]->create3(trellis->factory[j], l, sfl);
				f = zSFListMoveFirst(sfl);
				while(f != NULL){
					f->strand = '+';	
					f->state  = j;
					zSFListInsert(trellis->external[j], f);
					f = zSFListMoveNext(sfl);
				}
				zResetSFList(sfl);
			}
			l = trellis->dna->length - i - 1;
			if(trellis->rfactory[j] != NULL){
				zResetSFList(trellis->rexternal[j]);		
				trellis->rfactory[j]->create5(trellis->rfactory[j], l, sfl);
				f = zSFListMoveFirst(sfl);
				while(f != NULL){
					f->strand = '-';
					f->state  = j;
					tmp_start = f->start;
					f->start  = trellis->dna->length - f->end - 1;
					f->end    = trellis->dna->length - tmp_start - 1;
					zSFListInsert(trellis->rexternal[j], f);
					f = zSFListMoveNext(sfl);
				}
				zResetSFList(sfl);
			}
		}
			
		for (j = 0; j < hmm->states; j++) {
			if (NULL == trellis->cell[j]) continue;
			trellis->cell[j][i].score[0] = MIN_SCORE;
			trellis->cell[j][i].trace[0] = -1;
			trellis->cell[j][i].length[0] = 0;

			
			if ((hmm->state[j].type == EXPLICIT)
				|| (hmm->state[j].type == INTERNAL)){
				trellis->cell[j][i].submax = zMalloc(sizeof(zSFVec),
													 "zRunViterbiAndForward submax");
				zInitSFVec(trellis->cell[j][i].submax, nrange[j]);
				trellis->cell[j][i].submax->size = nrange[j];
				trellis->cell[j][i].submax->last =
					&trellis->cell[j][i].submax->elem[nrange[j]-1];
				
				for (u = 0; u < nrange[j]; u++){
					trellis->cell[j][i].submax->elem[u].score = MIN_SCORE;
					trellis->cell[j][i].submax->elem[u].intrinsic_score = MIN_SCORE;
					trellis->cell[j][i].submax->elem[u].from_state = -1;
					trellis->cell[j][i].submax->elem[u].start = i;
					trellis->cell[j][i].submax->elem[u].end = i;
				}
			}

			jumps = hmm->jmap[j];
			for(k = 0; k < jumps->size; k++) {
				zGetTransFunc(hmm->state[j].type) 
					(trellis, jumps->elem[k], j, i);
			}
			/*	if(!zOption("top")){ */
			if ((hmm->state[j].type == EXTERNAL)
			   	&& (trellis->cell[j][i].score[0] > MIN_SCORE)){ 
			   	zPushIVec(&trellis->extpos[j], i);
			}
			/*	}
				else{
				if(hmm->state[j].type == EXTERNAL){
					for(m=0;m<top_num;m++){
						if(trellis->cell[j][i].score[m] > MIN_SCORE){
							zPushIVec(&trellis->extpos[m][j]);
						}
					}
					}*/
		
			
			if ((hmm->state[j].type == EXPLICIT)

				|| (hmm->state[j].type == INTERNAL)){
				if (trellis->cell[j][i-1].submax->elem != NULL){
					zFreeSFVec(trellis->cell[j][i-1].submax);
					trellis->cell[j][i-1].submax->elem = NULL;
				}
				zFree(trellis->cell[j][i-1].submax);
				trellis->cell[j][i-1].submax = NULL;
			}
		}
	}
	
	zFreeSFList(sfl);
	zFree(sfl);

	/* Viterbi trace back */
	zTrace2("traceback\n");
	sfv = zCalloc(top_num,sizeof(zSFVec*),"zRunViterbiAndForward sfv");
	size = sizeof(zSFVec);
    for(g=0;g<top_num;g++){
		sfv[g]=zMalloc(size,"zMalloc sfv[i]");
		zInitSFVec(sfv[g],10);
	}
zTraceTrellis(trellis, -1, sfv, path_score);
    for(g=0;g<top_num;g++){
		zTranslateSFVec(sfv[g], -PADDING);
	}
	/*
 *path_score = trellis->cell[sfv->elem[0].state][dna->length-1-PADDING].score[0];
 */
	for (j = 0; j < hmm->states; j++){
		if ((hmm->state[j].type == EXPLICIT)
			|| (hmm->state[j].type == INTERNAL)){
			if (trellis->cell[j][dna->length-1].submax->elem != NULL){
				zFreeSFVec(trellis->cell[j][dna->length-1].submax);
				trellis->cell[j][dna->length-1].submax->elem = NULL;
			}
			zFree(trellis->cell[j][dna->length-1].submax);
			trellis->cell[j][dna->length-1].submax = NULL;
		}
	}
	zFree(nrange);
	return sfv;
}

void zShowTrellis(zTrellis *trellis){
	zHMM         *hmm = trellis->hmm;
	/*	zDNA         *dna = trellis->dna; */

	coor_t        i;        /* iterator for sequence */
	int           j, k;        /* iterator for internal states */

	j = 0;
	printf("===================== Show Trellis ===================================\n");
	/*for(i = PADDING; i< dna->length - PADDING; i ++) { */
	for(i = 48; i< 53; i ++) { 
		if(j==10){
			for(k= 0; k <  hmm->states; k++) {
				printf("%s\t", zStrIdx2Char(k));
			}
			printf("\n");		
			j = 0;
		}
		else{ j ++;}
		printf("%d\t", i-PADDING);
		for(k= 0; k < 4; k++) {

			printf("%3.0f\t",trellis->cell[k][i].score[0]);
		printf("***%3.0f\t",trellis->cell[k][i].score[1]);
			printf("###%d\t",trellis->cell[k][i].trace_top[0]);
			printf("***%d\t",trellis->cell[k][i].trace_top[1]);	
	}
		printf("\n");		
	}
	printf("======================= End Trellis =================================\n");

}
void zRunBackward(zTrellis *trellis) {
	zHMM*         hmm = trellis->hmm;
	coor_t        i;        /* iterator for sequence */
	coor_t        length = trellis->dna->length;
	int           j;        /* iterator for internal states */
	int           k;        /* iterator for previous states */
	
	/* Init */
	zAllocBackwardVars(trellis);
	for (j = 0; j < hmm->states; j++) {
		for (i = length - PADDING; i < length; i++) {
			trellis->backward[j][i] = 0;
		}
	}

	/* Induction */
	for (i = length - PADDING - 1; i > 1; i--) {
		for (j = 0; j < hmm->states; j++) {
			if (NULL == trellis->backward[j]) continue;
			for (k = 0; k < hmm->states; k++) {
				if (MIN_SCORE == zGetTransitionScore(trellis->hmm, j, k,trellis->tiso_group)) continue;
				zGetBackTransFunc(hmm->state[k].type)(trellis, j, k, i);
			}
		}
	} 
}

float  zForwardBackwardProb(zTrellis* trellis, zSfeature* f) {
	/*zHMM    *hmm = trellis->hmm;
	  zIVec   *states;
	  int      i, pre_state, state, post_state;
	  zPhase_t  phs;
	  score_t  parse_score, dscore, fscore, t1score, t2score, phscore;
	  score_t  pre_score, post_score;
	  score_t  result = MIN_SCORE;
	*/
	/* find the correct HMM state for the feature */
	/*	states = &hmm->simap[hmm->somap[f->name]];
	for (i = 0; i < states->size; i++) {
		state = states->elem[i];
		phs   = hmm->state[state].phase;
		if (zCompatibleJumpFromExon(trellis->dna, phs, f)) continue;
	}
	if (i == states->size) {
		zWarn("couldn't find the corresponding HMM state");
		return 0.0;
	}
	
	dscore = zScoreDuration(hmm->dmap[hmm->state[state].duration].score, 
							f->end - f->start + 1);
	fscore = f->score;
	*/
	/* sum over all paths that use this feature */
	/*	for (pre_state = 0; pre_state < hmm->states; pre_state++) {
		t1score = hmm->tmap[pre_state][post_state];
		pre_score = trellis->forward[pre_state][f->start-1];
		if (MIN_SCORE == t1score) continue;
		for (post_state = 0; post_state < hmm->states; post_state++) {
			t2score = hmm->tmap[pre_state][post_state];
			post_score = trellis->backward[post_state][f->end+1];
			phscore = zFhmm->phasepref
			Ack!! get the previous imp. from 8/15/02!
	*/
	trellis = NULL; /*compiler hush*/
	f = NULL;       /*compiler hush*/
	return 0.0;
}

score_t zGetIntergenicContinueScore(zTrellis* trellis) {
	int iso_group = trellis->tiso_group;

	if (iso_group < 0 || iso_group >= trellis->hmm->iso_transitions ) {
		zDie("zGetIntergenicContinueScore: iso_group %d not found\n",
			 iso_group);
	}
	return trellis->hmm->inter_continue[iso_group];
}

void zClearColumn(zTrellis* trellis, int pos) {
	int state;
	for (state=0; state < trellis->hmm->states; ++state) {
		trellis->cell[state][pos].length[0] = 0;
		trellis->cell[state][pos].score[0]  = MIN_SCORE;
		trellis->cell[state][pos].trace[0] = -1;
	}
}

score_t zGetPathScore(zTrellis* trellis, const zSFVec* sfv, zSFVec* result) {
	int             cur_state, prev_state;
	int             i, j;
	coor_t          begin, end, pos;
	zHMM           *hmm = trellis->hmm;
	zHMM_StateType  state_type;
	zSFVec         **external;/*, *exter;
							   zFeatureFactory *ff;*/
	score_t  *score;
	/* initialize trellis */
	end = PADDING-1;
	zAllocViterbiVars(trellis);
	for (pos = 0; pos <= end; ++pos) {
		zClearColumn(trellis, pos);
	}
	for (cur_state = 0; cur_state < hmm->states; cur_state++) {
		trellis->cell[cur_state][end].score[0] = zGetInitProb(trellis->hmm, cur_state, trellis->iiso_group);
		trellis->cell[cur_state][end].trace[0] = cur_state;
	}

	/* loop set up */
	if (sfv->size > 0) 
		cur_state = sfv->elem[0].state; /* for prev_state, below */

	/* trace through each feature */
	for (i = 0; i < sfv->size; ++i) {
		zSfeature      *f;
        external = zMalloc(sizeof(zSFVec), "external in zGetPathScore");
        zInitSFVec(external[0], 1);

		/* get vitals on this feature */
		f = &sfv->elem[i];
		prev_state = cur_state;
		cur_state  = f->state;
		begin = end + 1;

		/* run Viterbi if necessary */
		if (begin != f->start) {
			coor_t pos;
			for (pos = begin; pos < f->start; pos++) {
				zClearColumn(trellis, pos);
			}
			zRunPartialViterbiAndForward(trellis, begin, f->start-1);
			prev_state = -1;
		}
		end = f->end;
		
		/* transition from prev_state to cur_state */
		state_type = hmm->state[cur_state].type;
		if (INTERNAL == state_type || GINTERNAL == state_type ) { 
			/*Internal or Ginternal transitions go one base at a time */
			zClearColumn(trellis, begin);
			if (prev_state < 0) {
				zRunPartialViterbiAndForward(trellis, f->start, f->start);
			} else {
				zGetTransFunc(state_type)(trellis,prev_state,cur_state,begin);
			}
			for(pos = f->start+1; pos <= end; ++pos) {
				zClearColumn(trellis, pos);
				zGetTransFunc(state_type)(trellis, cur_state, cur_state, pos);
				trellis->cell[cur_state][pos-1].score[0] = MIN_SCORE;
				if (trellis->cell[cur_state][pos].score[0] < -10000) {
					printf("messed up on feature %d, at %u\n", i, pos);
				}
			}
		} else { /* Other states make jumps */
			for (pos = f->start; pos <= end; ++pos) {
				zClearColumn(trellis, pos);
			}
			/*EVAN I'm guessing that this comment breaks this so I will need to fix it later
			if (EXTERNAL == state_type) { / put only that feature in the vec/
				exter = ('-' == hmm->state[cur_state].strand) ?
					&trellis->rexternal[hmm->state[cur_state].model][end] :
					&trellis->external[hmm->state[cur_state].model][end];
				zFreeSFVec(exter);
				zInitSFVec(exter, 1);
				ff = ('-' == f->strand)
					? trellis->rfactory[hmm->state[cur_state].model]
					: trellis->factory[hmm->state[cur_state].model];
				if ('-' == f->strand) zAntiSfeature(f, trellis->dna->length);
				f->score = ff->scoref(ff, f);
				if ('-' == f->strand) zAntiSfeature(f, trellis->dna->length);
				zPushSFVec(exter, f);
			}
			*/
			if (prev_state < 0) {
				zRunPartialViterbiAndForward(trellis, end, end);
			} else {
				zGetTransFunc(state_type)(trellis, prev_state, cur_state, end);
			}
		}

		/* gather the results of the transitions */
		zTracePartialTrellis(trellis, begin, end, cur_state, external, score);
		if (NULL == external) {
			zWarn("Couldn't trace back from state %s, is the feature legal?",
				  zStrIdx2Char(hmm->state[cur_state].name));
			return MIN_SCORE;
		}
		for (j = external[0]->size; j > 0; j--) {
			zPushSFVec(result, &external[0]->elem[j-1]);
		}
		zFreeSFVec(external[0]);
		zFree(external[0]);
		external[0] = NULL;

		/* clean up so that it's impossible to jump back here */
		for (pos = begin-1; pos < end; pos++) {
			zClearColumn(trellis, pos);
		}

		/* and a sanity check */
		if (trellis->cell[cur_state][end].score[0] < -10000) {
			zWarn("Model prohibits ending in %s, at %d", 
				  zStrIdx2Char(hmm->state[f->state].name), 
				  end - trellis->padding + 1);
			return MIN_SCORE;
		}
	}

	/* run one last viterbi, if needed */
	if (end < trellis->dna->length - PADDING - 1) {
		coor_t pos;
		external = zMalloc(sizeof(zSFVec), "external in zGetPathScore");
		zInitSFVec(external[0], 1); 

		begin = end + 1;
		end = trellis->dna->length - PADDING - 1;
		for (pos = begin; pos < end; pos++) {
			zClearColumn(trellis, pos);
		}
		zRunPartialViterbiAndForward(trellis, begin, end);
		zTracePartialTrellis(trellis, begin, end, -1, external,score);
		for (j = external[0]->size; j > 0; j--) {
			zPushSFVec(result, &external[0]->elem[j-1]);
		}
		zFreeSFVec(external[0]);
		zFree(external[0]);
		external[0] = NULL;
	}
	assert(end == trellis->dna->length - PADDING - 1);

	/* report end score */
	return trellis->cell[result->last->state][end].score[0];
}

frame_t zFragFrame2Char(frame_t lfrag, frame_t rfrag, frame_t frame)
{
	/* Encodes lfrag, rfrag and frame values in a character      */
	/* in bitwise manner, as shown below. All values range       */
	/* between -1 and 2, hence 2 bits are enough to encode each. */

	/* --------------------------------------------------------- */
	/* bit index          :   7 6     5 4     3 2     1 0        */
	/*                                                           */
	/* character (1 byte) :   0 0 --- 0 1 --- 1 0 --- 0 0        */
	/*                        ^ ^     ^ ^     ^ ^     ^ ^        */
	/*                        | |     | |     | |     | |        */
	/* encoded info       :  unused  frame   rfrag   lfrag       */
	/*                                                           */
	/* values             :            1       2       0         */
	/* --------------------------------------------------------- */

	frame_t result = 0; /* frame_t is typedef'ed as char */

	if ((lfrag > 2) || (rfrag > 2) || (frame > 2))
		{
			fprintf(stderr, "lfrag, rfrag or frame is out of bounds: %d %d %d\n",
					lfrag, rfrag, frame);
		}
	else
		{
			/* value of -1 is allowed for promoters */
			/* and PolyAs and is encoded here as 3  */
			if (lfrag < 0) lfrag = 3;
			if (rfrag < 0) rfrag = 3;
			if (frame < 0) frame = 3;

			result |= lfrag;        /* set bits 0 and 1 using lfrag */

			result |= (rfrag << 2); /* set bits 2 and 3 using rfrag */
			/* shifted left by 2 bits       */

			result |= (frame << 4); /* set bits 4 and 5 using frame */
			/* shifted left by 4 bits       */
		}

	return(result);
}

void zChar2FragFrame(frame_t frag_frame, zSfeature *f)
{
	/* Reads lfrag, rfrag and frame values encoded in character frag_frame */
	/* and assigns those to the corresponding members of zSfeature object. */
	/* The bytes are read two at a time by means of bitwise AND operation  */
	/* with a number that has 1's in the bits to be read and 0's elsewhere */

  /* ------------------------------------------------------------------- */
  /* bit index          :   7 6     5 4     3 2     1 0                  */
  /*                                                                     */
  /* character (1 byte) :   0 0 --- 0 1 --- 1 0 --- 0 0                  */
  /*                                                & &                  */
  /*     & 3            :                           1 1 = 0 0 (lfrag=0)  */
  /*                                        & &                          */
  /*     & (3 << 2)     :                   1 1 >> 2    = 1 0 (rfrag=2)  */
  /*                                & &                                  */
  /*     & (3 << 4)     :           1 1 >> 4            = 0 1 (frame=1)  */
  /* ------------------------------------------------------------------- */

	f->lfrag = frag_frame & 3;               /* set lfrag from bits 0 and 1 */
	f->rfrag = (frag_frame & (3 << 2)) >> 2; /* set rfrag from bits 2 and 3 */
	f->frame = (frag_frame & (3 << 4)) >> 4; /* set frame from bits 4 and 5 */

  
	/* value of -1 is allowed for promoters */
	/* and PolyAs and is encoded here as 3  */
	if (f->lfrag == 3) f->lfrag = -1;
	if (f->rfrag == 3) f->rfrag = -1;
	if (f->frame == 3) f->frame = -1;
}


/**********************************************************\
  Pin Data Search Code                                   
	EVAN
\**********************************************************/

/* EVAN removed checking code because it is old style, need to fix it and test feature creation 
static void check_for_feature(zSfeature *f,zSFList *sfl,int erase){
	zSfeature *f2;
	f2 = zSFListMoveFirst(sfl);
	while(f2 != NULL){
		if(f->start == f->start && f->end == f2->end && f->score == f2->score && f->lfrag == f2->lfrag && f->rfrag == f2->rfrag && f->frame == f2->frame){
			if(erase != 0){
				zSFListRemoveCurrent(sfl);
			}
			return;
		}
		f2 = zSFListMoveNext(sfl);
	}
	if(erase){
		fprintf(stderr,"Missed %d (%d->%d) %c\n",f->name,f->start,f->end,f->strand);
	}
	else{
		fprintf(stderr,"Created bad %d (%d->%d) %c\n",f->name,f->start,f->end,f->strand);
	}
}

static void check_for_all_features(zTrellis *trellis,zSFList *sfl, coor_t pos, int fnum, strand_t strand){
	coor_t i,end;
	zSFList* sfl2;
	zSfeature *f2;
	
	end = pos + 50000;
	if(end > trellis->dna->length) end = trellis->dna->length -10;
	for(i = pos; i < end;i++){
		if(strand == '+'){
			sfl2 = zGetTrellisExternalValue(trellis,i,fnum);
		}
		else{
			sfl2 = zGetTrellisRExternalValue(trellis,i,fnum);
		}
		f2 = zSFListMoveFirst(sfl2);
		while(f2 != NULL){
			if(f2->start < pos){
				check_for_feature(f2,sfl,1);
			}
			f2 = zSFListMoveNext(sfl2);
		}
	}
}

static void check_inside_create(zTrellis *trellis){
	int j;
	coor_t i;
	zSFList sfl;
	zSfeature *f;

	zInitSFList(&sfl);

	zComputeBackLinks2(trellis);
	dprintf(0,"External vals ready\n");
	for(i = 10;i <= trellis->dna->length-10; i++){
		if(i % 100 == 0){
			dprintf(9,"Checked up to %d\n",i);
		}
		for(j = 0;j < trellis->hmm->feature_count;j++){
			if(trellis->factory[j] != NULL){
				trellis->factory[j]->icreate(trellis->factory[j],i,&sfl);
				f = zSFListMoveFirst(&sfl);
				while(f != NULL){
					f->strand = '+';
					check_for_feature(f,zGetTrellisExternalValue(trellis,f->end,j),0);
					f = zSFListMoveNext(&sfl);
				}
				check_for_all_features(trellis,&sfl,i,j,'+');
				zResetSFList(&sfl);
			}
			if(trellis->rfactory[j] != NULL){
				trellis->factory[j]->icreate(trellis->rfactory[j],i,&sfl);
				f = zSFListMoveFirst(&sfl);
				while(f != NULL){
					f->strand = '-';
					check_for_feature(f,zGetTrellisRExternalValue(trellis,f->end,j),0);
					f = zSFListMoveNext(&sfl);
				}
				check_for_all_features(trellis,&sfl,i,j,'-');
				zResetSFList(&sfl);
			}
		}
	}
}

static void check_link_dist(zTrellis *trellis) {
	coor_t     i;
	int        j;
	coor_t     tmp_start;
	zSFList   *sfl;
	zSfeature *f;

	int idx;
	coor_t split_size = 5000;
	int    split_count = trellis->dna->length/split_size + 1;
	int*   count = zMalloc(sizeof(int)*split_count,"");
	for(j = 0;j < split_count;j++){
		count[j] = 0;
	}

	if (NULL != trellis->external) return;

	trellis->external = zMalloc(sizeof(zExternal),
															"zComputeBackLinks trellis->external");
	zInitExternal(trellis->external,trellis->factory,1,
								trellis->hmm->feature_count);
	trellis->rexternal = zMalloc(sizeof(zExternal),
															 "zComputeBackLinks trellis->rexternal");
	zInitExternal(trellis->rexternal,trellis->rfactory,1,
								trellis->hmm->feature_count);
	sfl = zMalloc(sizeof(zSFList), "zComputeBackLinks sfl");
	zInitSFList(sfl);
  
	for (j = 0; j < trellis->hmm->feature_count; j++) {
		if(trellis->factory[j] != NULL){
			for (i = 0; i < trellis->dna->length; i++) {
				trellis->factory[j]->create(trellis->factory[j], i, sfl);
				f = zSFListMoveFirst(sfl);
				while(f != NULL){
					f->strand = '+';	

					idx = (int)f->start/split_size;
					while(idx*split_size <= f->end){
						if(idx*split_size >= f->start){
							count[idx]++;
						}
						idx++;
					}
					f = zSFListMoveNext(sfl);
				}
				zResetSFList(sfl);
			}
		}
		if(trellis->rfactory[j] != NULL){
			for (i = 0; i < trellis->dna->length; i++) {		
				trellis->rfactory[j]->create(trellis->rfactory[j], i, sfl);
				f = zSFListMoveFirst(sfl);
				while(f != NULL){
					f->strand = '-';
					tmp_start = f->start;
					f->start  = trellis->dna->length - f->end - 1;
					f->end    = trellis->dna->length - tmp_start - 1;

					idx = (int)f->start/split_size;
					while(idx*split_size <= f->end){
						if(idx*split_size >= f->start){
							count[idx]++;
						}
						idx++;
					}
					f = zSFListMoveNext(sfl);
				}
				zResetSFList(sfl);
			}
		}
	}
  
	zFreeSFList(sfl);
	zFree(sfl);

	for(j = 0;j < split_count;j++){
		dprintf(-1000,"pos\t%d\t%d\n",j*split_size,count[j]);
	}
}
*/

struct zPinPathData{
	score_t score;
	int     pinstate;
};
typedef struct zPinPathData zPinPathData;

static void zInitPinPathData(zPinPathData* p){
	p->score = MIN_SCORE;
	p->pinstate = -1;
}

struct zPinArrayData{
	int    index;
	int    state;
	coor_t pos;
	/*the following are needed only for 5' long distance links*/
	int link_state;
};
typedef struct zPinArrayData zPinArrayData;

static void* zCreatePinArrayData(){
	zPinArrayData* pad = zMalloc(sizeof(zPinArrayData),"zCreateInt pad");
	pad->index = -1;
	pad->state = -1;
	pad->pos = (coor_t)-1;
	return (void*)pad;
}

static void zFreePinArrayData(void* v){
	zFree(v);
}

static void zResetPinArrayData(void* v){
	zPinArrayData* pad = (zPinArrayData*)v;
	pad->index = -1;
	pad->state = -1;
	pad->pos = (coor_t)-1;
}

static void* zCreateInt(){
	int* i = zMalloc(sizeof(int),"zCreateInt i");
	return i;
}

static void zFreeInt(void* v){
	zFree(v);
}

static void zResetInt(void* v){
	int* i = (int*)v;
	*i = 0;
}

static void zInitIList(zList* l){
	zInitList(l,zCreateInt,zFreeInt,zResetInt);
}

static void zIListAddLast(zList* l, int i){
	int *j;
	j = zListAddLast(l);
	*j = i;
}

static void zIListAddFirst(zList* l, int i){
	int *j;
	j = zListAddFirst(l);
	*j = i;
}

static int zIListPopFirst(zList* l){
	int *j;
	j = zListMoveFirst(l);
	zListRemoveFirst(l);
	if(j == NULL){
		zDie("Poped empty zIList");
	}
	return *j;
}

struct zPinSearch{
	zTrellis*          trellis;          /* zTrellis object */
	zPinPathData**     data;             /* matrix holding all path data */
	int                states;           /* number of states */
	int                size;             /* data matrix size */
	zList              live_to_list;     /* list of live path start nodes
											maps path start to zPinSearch.data row
											sorted by decreasing seq pos */
	zList              live_from_list;   /* list of live path end nodes
											maps path end to zPinSearch.data column
											sorted by increasing seq pos */
	zList              free_to_list;     /* list of unused rows in zPinSearch.data */
	zList              free_from_list;   /* list of unused columns in zPinSearch.data */
	zPinPathData*      tmp;               /* temp buffer to hold data when stepping 
											 forward of backward by one pos */
	coor_t             pin_pos;          /* the pos at which we are searching for a pin */
	coor_t             start_pos;        /* the current maximum path start pos */
	coor_t             end_pos;          /* the current minimum path end pos 
											at any time we know all path which can get us 
											from start_pos to end_pos, including any that 
											include a long distance link across either 
											start_pos or end_pos */
	int                pin_state;        /* the pin state at pin_pos, this is -1 until we
											find the correct value */
	zSFList            sfl;              /* a temp zSfeature list to avoid repeated 
											allocation and initialization */ 
};
typedef struct zPinSearch zPinSearch;

static int INIT_PIN_SEARCH_SIZE = 2500;
static const coor_t PIN_SEARCH_STEPS = 5000;

static void zStepPinSearchBack(zPinSearch*);
static void zStepPinSearchForward(zPinSearch*);
static void zPinSearchCheckForPin(zPinSearch*);
static void zStepPinArrayForward(zPinSearch*, zPinArrayData*);
static void zStepPinArrayBack(zPinSearch*, zPinArrayData*);
static void zPinSearchCreateLiveNodes3(zPinSearch*);
static void zPinSearchCreateLiveNodes5(zPinSearch*);
static void zPinSearchAddFeature3(zPinSearch*, int, strand_t);
static void zPinSearchAddFeature5(zPinSearch*, int, strand_t);
static void zPinSearchAddFeatureInit(zPinSearch*, int, strand_t);
static zPinArrayData* zPinSearchGetToPad(zPinSearch*, coor_t, int);
static zPinArrayData* zPinSearchGetFromPad(zPinSearch*, coor_t, int);

/* EVAN currently ignoring the long distance links which cross the pin_pos */
/* EVAN stupidly adding multiple list nodes for the same start/end state-pos combo */

static void zInitPinSearch(zPinSearch* ps, zTrellis* trellis, coor_t pin_pos){
	int i,j;
	zPinArrayData* pad;
	ps->trellis = trellis;
	ps->states = trellis->hmm->states;
	ps->pin_pos = pin_pos;
	ps->start_pos = pin_pos;
	ps->end_pos = pin_pos;
	ps->pin_state = -1;

	/*EVAN no current mechanism for increasing ps->size.  If it isn't big enough
		program will exit with error via zDie()*/
	ps->size = INIT_PIN_SEARCH_SIZE;
	ps->data = zMalloc(sizeof(zPinPathData*)*ps->size,"zFindTrellisPin data");
	ps->tmp = zMalloc(sizeof(zPinPathData)*ps->size,"zFindTrellisPin tmp");
	for(i = 0;i < ps->size;i++){
		zInitPinPathData(&ps->tmp[i]);
		ps->data[i] = zMalloc(sizeof(zPinPathData)*ps->size,"zFindTrellisPin data[i]");
		for(j = 0;j < ps->size;j++){
			zInitPinPathData(&ps->data[i][j]);
		}
	}

	zInitList(&ps->live_to_list,zCreatePinArrayData,zFreePinArrayData,zResetPinArrayData);
	zInitList(&ps->live_from_list,zCreatePinArrayData,zFreePinArrayData,
						zResetPinArrayData);
	zInitIList(&ps->free_to_list);
	zInitIList(&ps->free_from_list);
	
	/* Add row/col to live list for each state, i, at pos i+1 (adjusted for tmp at 0) */
	for(i = 0;i < ps->states;i++){
		pad = zListAddLast(&ps->live_to_list);
		pad->index = i;
		pad->pos = pin_pos;
		pad->state = i;
		pad = zListAddLast(&ps->live_from_list);
		pad->index = i;
		pad->pos = pin_pos;
		pad->state = i;
		
		if(trellis->hmm->state[i].type == INTERNAL ||
		   trellis->hmm->state[i].type == GINTERNAL ){
			ps->data[i][i].pinstate = i;
			ps->data[i][i].score = 0;
		}
		else if(trellis->hmm->state[i].type == EXTERNAL){
			ps->data[i][i].pinstate = -1;
			ps->data[i][i].score = MIN_SCORE;
		}
		else{
			zDie("Pin Search does not support states other than INTERNAL and EXTERNAL");
		}
	}
	/* Add the rest of the rows/cols to the free list */ 
	for(i = ps->states;i < ps->size-1;i++){
		zIListAddLast(&ps->free_to_list,i);
		zIListAddLast(&ps->free_from_list,i);
	}

	zInitSFList(&ps->sfl);

	/* Add long distance links crossing pin_pos */
	for (i = 0; i < ps->trellis->hmm->feature_count; i++) {
		j = ps->pin_pos;
		if(trellis->factory[i] != NULL){
			/*dprintf(0,"Adding Feature %d at pos %d +\n",i,j);*/
			trellis->factory[i]->createi(trellis->factory[i], j, &ps->sfl);
			zPinSearchAddFeatureInit(ps,i,'+');
			zResetSFList(&ps->sfl);
		}
		j = trellis->dna->length - 1 - ps->pin_pos;
		if(trellis->rfactory[i] != NULL){
			/*dprintf(0,"Adding Feature %d at pos %d -\n",i,j);*/
			trellis->rfactory[i]->createi(trellis->rfactory[i], j, &ps->sfl);
			zPinSearchAddFeatureInit(ps,i,'-');
			zResetSFList(&ps->sfl);
		}		
	}	
	fprintf(stderr,"Done adding Features\n");

}

static void zFreePinSearch(zPinSearch* ps){
	int i;
	zFreeSFList(&ps->sfl);
	zFreeList(&ps->live_to_list);
	zFreeList(&ps->live_from_list);
	zFreeList(&ps->free_to_list);
	zFreeList(&ps->free_from_list);
	zFree(ps->tmp);
	ps->tmp = NULL;
	for(i = 0;i < ps->size;i++){
		zFree(ps->data[i]);
		ps->data[i] = NULL;
	}
	zFree(ps->data);
	ps->data = NULL;
}

static void zPinSearchResetTo(zPinSearch* ps, int i){
	zPinArrayData* pad;
	if(i >= ps->size){
		return;
	}
	pad = zListMoveFirst(&ps->live_from_list);
	while(pad != NULL){
		ps->data[pad->index][i].score = MIN_SCORE;
		ps->data[pad->index][i].pinstate = -1;
		pad = zListMoveNext(&ps->live_from_list);
	}
}

static void zPinSearchResetFrom(zPinSearch* ps, int i){
	zPinArrayData* pad;
	if(i >= ps->size){
		return;
	}
	pad = zListMoveFirst(&ps->live_to_list);
	while(pad != NULL){
		ps->data[i][pad->index].score = MIN_SCORE;
		ps->data[i][pad->index].pinstate = -1;
		pad = zListMoveNext(&ps->live_to_list);
	}
}

static void zPinSearchResetTmp(zPinSearch* ps){
	int i;
	for(i = 0;i < ps->states;i++){
		ps->tmp[i].score = MIN_SCORE;
		ps->tmp[i].pinstate = -1;
	}
}

static zPinArrayData* zPinSearchGetToPad(zPinSearch* ps, coor_t pos, int state){
	zPinArrayData* to_pad;

	to_pad = zListGetCurrent(&ps->live_to_list);
	if(to_pad == NULL){
		/* start at end of list (highest coor)*/
		to_pad = zListMoveLast(&ps->live_to_list);
	}
	/* if the list isn't empty */
	if(to_pad != NULL){
		if(to_pad->pos > pos){
			/* move back until we are at a node <= pos */
			while(to_pad != NULL && to_pad->pos > pos){
				to_pad = zListMovePrev(&ps->live_to_list);
			}
			/* while the current node == pos, see if it is the correct state */
			while(to_pad != NULL && to_pad->pos == pos){
				if(to_pad->state == state){
					return to_pad;
				}
				to_pad = zListMovePrev(&ps->live_to_list);
			}
			/* didn't find the node in the list so add a new one */
			to_pad = zListAddNext(&ps->live_to_list);
		}
		else{
			/* move forward until we are at a node >= pos */
			while(to_pad != NULL && to_pad->pos < pos){
				to_pad = zListMoveNext(&ps->live_to_list);
			}
			/* while the current node == pos, see if it is the correct state */
			while(to_pad != NULL && to_pad->pos == pos){
				if(to_pad->state == state){
					return to_pad;
				}
				to_pad = zListMoveNext(&ps->live_to_list);
			}
			/* didn't find the node in the list so add a new one */
			to_pad = zListAddPrev(&ps->live_to_list);
		}
	}
	/* fill in the new nodes info */
	to_pad->index = zIListPopFirst(&ps->free_to_list);
	to_pad->state = state;
	to_pad->pos = pos;
	/*dprintf(0,"ADD TOPAD %d %d\n",pos,state);EVAN*/
	return to_pad;
}

static zPinArrayData* zPinSearchGetFromPad(zPinSearch* ps, coor_t pos, int state){
	zPinArrayData* from_pad;

	from_pad = zListGetCurrent(&ps->live_from_list);
	if(from_pad == NULL){
		/* start at end of list (lowest coor)*/
		from_pad = zListMoveLast(&ps->live_from_list);
	}
	/* if the list isn't empty */
	if(from_pad != NULL){
		if(from_pad->pos < pos){
			/* move back until we are at a node >= pos */
			while(from_pad != NULL && from_pad->pos < pos){
				from_pad = zListMovePrev(&ps->live_from_list);
			}
			/* while the current node == pos, see if it is the correct state */
			while(from_pad != NULL && from_pad->pos == pos){
				if(from_pad->state == state){
					return from_pad;
				}
				from_pad = zListMovePrev(&ps->live_from_list);
			}
			/* didn't find the node in the list so add a new one */
			from_pad = zListAddNext(&ps->live_from_list);
		}
		else{
			/* move forward until we are at a node <= pos */
			while(from_pad != NULL && from_pad->pos > pos){
				from_pad = zListMoveNext(&ps->live_from_list);
			}
			/* while the current node == pos, see if it is the correct state */
			while(from_pad != NULL && from_pad->pos == pos){
				if(from_pad->state == state){
					return from_pad;
				}
				from_pad = zListMoveNext(&ps->live_from_list);
			}
			/* didn't find the node in the list so add a new one */
			from_pad = zListAddPrev(&ps->live_from_list);
		}
	}
	/* fill in the new nodes info */
	from_pad->index = zIListPopFirst(&ps->free_from_list);
	from_pad->state = state;
	from_pad->pos = pos;
	/*dprintf(0,"ADD FROMPAD %d %d\n",pos,state);EVAN*/
	return from_pad;
}

void zCheckForwardSearch(zPinSearch* ps){
	coor_t i;
	int keep_state = 0;
	zPinArrayData *to_pad,*from_pad;
	coor_t lead_in = 100000;

	ps->start_pos -= lead_in;
	ps->end_pos -= lead_in;
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL){
		to_pad->pos -= lead_in;
		to_pad = zListMoveNext(&ps->live_to_list);
	}
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL){
		from_pad->pos -= lead_in;

		if(from_pad->state == 9999){/*keep_state){*/	
			/*remove long distance link node */
			/*dprintf(0,"\tremoving %d\n",from_pad->state);*/
			zPinSearchResetFrom(ps,from_pad->index);
			zIListAddFirst(&ps->free_from_list,from_pad->index);
			zListRemoveCurrent(&ps->live_from_list);
		}

		from_pad = zListMoveNext(&ps->live_from_list);
	}

	if(ps->start_pos > ps->trellis->dna->length){
		zDie("I suck\n");
	}

	for(i = 0;i < lead_in;i++){
		zStepPinSearchForward(ps);
	}
	fprintf(stderr,"Search forward to %d (%d->%d)\n",ps->end_pos,ps->live_from_list.size,ps->live_to_list.size);
	zPinSearchCheckForPin(ps);

	fprintf(stderr,"Setting pin at pos %d\n",ps->end_pos);
	fprintf(stderr,"trimming to only state %d at start (%d)\n",keep_state,ps->start_pos);
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL){
		to_pad = zListMoveFirst(&ps->live_to_list);
		while(to_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score != MIN_SCORE){
				/*dprintf(0,"\t\t(%d->%d) (%d->%d) pin = %d->%d %f\n",from_pad->state,to_pad->state,from_pad->pos,to_pad->pos,ps->data[from_pad->index][to_pad->index].pinstate,to_pad->state,ps->data[from_pad->index][to_pad->index].score);			*/
				ps->data[from_pad->index][to_pad->index].pinstate = to_pad->state;			
			}
			to_pad = zListMoveNext(&ps->live_to_list);
		}
		from_pad = zListMoveNext(&ps->live_from_list);
	}
	fprintf(stderr,"Search forward to %d (%d->%d)\n",ps->end_pos,ps->live_from_list.size,ps->live_to_list.size);
	zPinSearchCheckForPin(ps);
	while(ps->pin_state < 0){
		for(i = 0;i < PIN_SEARCH_STEPS;i++){
			zStepPinSearchForward(ps);
		}
		fprintf(stderr,"Search forward to %d (%d->%d)\n",ps->end_pos,ps->live_from_list.size,ps->live_to_list.size);
		zPinSearchCheckForPin(ps);
	}
}


void zCheckBackSearch(zPinSearch* ps){
	coor_t i;
	int keep_state = 0;
	zPinArrayData *to_pad,*from_pad;
	coor_t lead_in = 100000;

	ps->start_pos += lead_in;
	ps->end_pos += lead_in;
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL){
		to_pad->pos += lead_in;
		if(to_pad->state == 9999){/*keep_state){*/
			zPinSearchResetTo(ps,to_pad->index);
			zIListAddFirst(&ps->free_to_list,to_pad->index);
			zListRemoveCurrent(&ps->live_to_list);
		}
		to_pad = zListMoveNext(&ps->live_to_list);
	}
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL){
		from_pad->pos += lead_in;

		from_pad = zListMoveNext(&ps->live_from_list);
	}

	if(ps->end_pos > ps->trellis->dna->length){
		zDie("I suck\n");
	}
	
	for(i = 0;i < lead_in;i++){
		zStepPinSearchBack(ps);
	}
	fprintf(stderr,"Search back to %d (%d->%d)\n",ps->start_pos,ps->live_from_list.size,ps->live_to_list.size);
	zPinSearchCheckForPin(ps);

	fprintf(stderr,"Setting pin at pos %d\n",ps->start_pos);
	fprintf(stderr,"trimming to only state %d at start (%d)\n",keep_state,ps->end_pos);
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL){
		from_pad = zListMoveFirst(&ps->live_from_list);
		while(from_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score != MIN_SCORE){
				/*fprintf(stderr,"\t\t(%d->%d) (%d->%d) pin = %d->%d %f\n",from_pad->state,to_pad->state,from_pad->pos,to_pad->pos,ps->data[from_pad->index][to_pad->index].pinstate,to_pad->state,ps->data[from_pad->index][to_pad->index].score);*/
				if(from_pad->pos < ps->start_pos){
					ps->data[from_pad->index][to_pad->index].pinstate = from_pad->link_state;			
					
				}
				else{
					ps->data[from_pad->index][to_pad->index].pinstate = from_pad->state;			
				}
			}
			from_pad = zListMoveNext(&ps->live_from_list);
		}
		to_pad = zListMoveNext(&ps->live_to_list);
	}
	fprintf(stderr,"Search back to %d (%d->%d)\n",ps->start_pos,ps->live_from_list.size,ps->live_to_list.size);
	zPinSearchCheckForPin(ps);
	while(ps->pin_state < 0){
		for(i = 0;i < PIN_SEARCH_STEPS;i++){
			zStepPinSearchBack(ps);
		}
		fprintf(stderr,"Search back to %d (%d->%d)\n",ps->start_pos,ps->live_from_list.size,ps->live_to_list.size);
		zPinSearchCheckForPin(ps);
	}
}



int zFindTrellisPin(zTrellis* trellis, coor_t pin_pos){
	coor_t i;
	zDNA*        dna;
	zPinSearch*  ps;
	int          pin_state = -1;
	
	pin_pos += PADDING - 1;

	dna = trellis->dna;
	if(pin_pos >= dna->length){
		zDie("zFindTrellisPin pin_pos must be less than the sequence length");
	}
	if(pin_pos == 0){
		zDie("zFindTrellisPin pin_pos must be greater than 0");
	}
	
	/* EVAN testing code
	if(pin_pos == 0){
		check_link_dist(trellis);
		check_inside_create(trellis);
		return -1;
	}
	*/

	/* init pin data structure */
	ps = zMalloc(sizeof(zPinSearch),"zFindTrellisPin ps");
	zInitPinSearch(ps,trellis,pin_pos);
	
	/* Until we have found the pin */
	zPinSearchCheckForPin(ps);
	/*EVAN*/
	/*zCheckBackSearch(ps);*/
	while(ps->pin_state < 0){
		for(i = 0;i < PIN_SEARCH_STEPS;i++){
			if(ps->end_pos >= trellis->dna->length - 1){
				zDie("Pin Search past end of sequence");
			}
			zStepPinSearchForward(ps);
		}
		fprintf(stderr,"Search forward to %d (%d->%d)\n",ps->end_pos,ps->live_from_list.size,ps->live_to_list.size);
		zPinSearchCheckForPin(ps);
		for(i = 0;i < PIN_SEARCH_STEPS;i++){
			if(ps->start_pos == 0){
				zDie("Pin Search past start of sequence");
			}
			zStepPinSearchBack(ps);
		}
		fprintf(stderr,"Search back to %d (%d->%d)\n",ps->start_pos,ps->live_from_list.size,ps->live_to_list.size);
		zPinSearchCheckForPin(ps);
	}
	pin_state = ps->pin_state;
	
	/* free pin data */ 
	zFreePinSearch(ps);
	zFree(ps);

	/* free trellis */
	zFreeTrellis(trellis);
	zFree(trellis);

	/* Return the pin state */
	/* EVAN what if ps->pin_pos != pin_pos */
	return pin_state;
}


static void zPinSearchCheckForPin(zPinSearch* ps){
	int i;
	int paths[100];
	int max_path_to_state[100];
	int min_path_to_state[100];
	int max_path_from_state[100];
	int min_path_from_state[100];
	score_t max_path_score[100];
	score_t min_path_score[100];
	int state = -1;
	coor_t last_from = -1;
	coor_t last_to = -1;
	int live_paths = 0;
	zPinArrayData* to_pad = zListMoveFirst(&ps->live_from_list);
	zPinArrayData* from_pad = zListMoveFirst(&ps->live_from_list);

	for(i = 0;i <= ps->states;i++){
		paths[i] = 0;
		max_path_score[i] = MIN_SCORE;
		min_path_score[i] = MAX_SCORE;
		max_path_to_state[i] = -1;
		min_path_to_state[i] = -1;
		max_path_from_state[i] = -1;
		min_path_from_state[i] = -1;
	}

	if(to_pad != NULL && from_pad != NULL){
		state = ps->data[from_pad->index][to_pad->index].pinstate;
		last_from = from_pad->pos;
		last_to = to_pad->pos;
	}
	else{
		zDie("Empty live_to_list or live_from_list");
	}

	if(from_pad->pos != ps->start_pos){
		zDie("Live from list not starting at start pos\n");
	}
	
	while(from_pad != NULL){
		to_pad = zListMoveFirst(&ps->live_to_list);
		if(to_pad->pos != ps->end_pos){
			zDie("Live to list not starting at end pos\n");
		}
		last_to = to_pad->pos;
		while(to_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score > MIN_SCORE){
				i = ps->data[from_pad->index][to_pad->index].pinstate;
				live_paths++;
				paths[i]++;
				/*dprintf(2,"(%d,%d)->(%d,%d) = %d (%.4f)\tp[%d]=%d\n",from_pad->pos,from_pad->state,to_pad->pos,to_pad->state, ps->data[from_pad->index][to_pad->index].pinstate,ps->data[from_pad->index][to_pad->index].score,ps->data[from_pad->index][to_pad->index].pinstate,paths[ps->data[from_pad->index][to_pad->index].pinstate]);*/
				if(state != i){
					state = -1;
					ps->pin_state = -1;
					/*return;*/
				}
				if(max_path_score[i] < ps->data[from_pad->index][to_pad->index].score){
					max_path_score[i] = ps->data[from_pad->index][to_pad->index].score;
					max_path_to_state[i] = to_pad->state;
					max_path_from_state[i] = from_pad->state;
				}
				if(min_path_score[i] > ps->data[from_pad->index][to_pad->index].score){
					min_path_score[i] = ps->data[from_pad->index][to_pad->index].score;
					min_path_to_state[i] = to_pad->state;
					min_path_from_state[i] = from_pad->state;
				}
			}
			if(last_to > to_pad->pos){
				zDie("Out of order to");
			}
			last_to = to_pad->pos;
			to_pad = zListMoveNext(&ps->live_to_list);
		}
		if(last_from < from_pad->pos){
			zDie("Out of order from");
		}
		last_from = from_pad->pos;
		from_pad = zListMoveNext(&ps->live_from_list);
	}
	ps->pin_state = state;
	
	fprintf(stderr,"\tCHECK FOR_PIN: %d live paths\n",live_paths);
	for(i = 0;i <= ps->states;i++){
		if(paths[i] > 0){
			fprintf(stderr,"\tstate %d - %d\t%f (%d->%d)\t%f (%d->%d)\n",i,paths[i],min_path_score[i],min_path_from_state[i],min_path_to_state[i],max_path_score[i],max_path_from_state[i],max_path_to_state[i]);
		}
	}
}

static void zCheckPinSearchListOrder(zPinSearch* ps){
	zPinArrayData*  to_pad;
	zPinArrayData*  from_pad;
	coor_t last;
	to_pad = zListMoveFirst(&ps->live_to_list);
	last = to_pad->pos;
	while(to_pad != NULL){
		if(last > to_pad->pos){
			to_pad = zListMoveFirst(&ps->live_to_list);
			last = to_pad->pos;
			while(to_pad != NULL){
				if(last > to_pad->pos){
					fprintf(stderr,"\t**");
				}
				fprintf(stderr,"\t%d\t%d\n",to_pad->pos,to_pad->state);
				last = to_pad->pos;
				to_pad = zListMoveNext(&ps->live_to_list);
			}
			fprintf(stderr,"LAST\n");
			to_pad = zListMoveLast(&ps->live_to_list);
			last = to_pad->pos;
			while(to_pad != NULL){
				if(last < to_pad->pos){
					fprintf(stderr,"\t**");
				}
				fprintf(stderr,"\t%d\t%d\n",to_pad->pos,to_pad->state);
				last = to_pad->pos;
				to_pad = zListMovePrev(&ps->live_to_list);
			}

			
			zDie("WTF to\n");
		}
		last = to_pad->pos;
		to_pad = zListMoveNext(&ps->live_to_list);
	}
	from_pad = zListMoveFirst(&ps->live_from_list);
	last = from_pad->pos;
	while(from_pad != NULL){
		if(last < from_pad->pos){
			from_pad = zListMoveFirst(&ps->live_from_list);
			last = from_pad->pos;
			while(from_pad != NULL){
				if(last < from_pad->pos){
					fprintf(stderr,"\t**");
				}
				fprintf(stderr,"\t%d\t%d\n",from_pad->pos,from_pad->state);
				last = from_pad->pos;
				from_pad = zListMoveNext(&ps->live_from_list);
			}
			fprintf(stderr,"LAST\n");
			from_pad = zListMoveLast(&ps->live_from_list);
			last = from_pad->pos;
			while(from_pad != NULL){
				if(last > from_pad->pos){
					fprintf(stderr,"\t**");
				}
				fprintf(stderr,"\t%d\t%d\n",from_pad->pos,from_pad->state);
				last = from_pad->pos;
				from_pad = zListMovePrev(&ps->live_from_list);
			}

			
			zDie("WTF from\n");
		}
		last = from_pad->pos;
		from_pad = zListMoveNext(&ps->live_from_list);
	}
}

static void zStepPinSearchForward(zPinSearch* ps){
	zPinArrayData*  to_pad;
	zPinArrayData*  from_pad;

	/* Move end forward */
	ps->end_pos++;

	/* create new 3' live nodes (links starting at end_pos)*/
	zPinSearchCreateLiveNodes3(ps);

	/* foreach live from-node update path to each live to-node at time end_pos */
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL){
		zStepPinArrayForward(ps,from_pad);
		from_pad = zListMoveNext(&ps->live_from_list);
	}

	/* update live_to_list states pad->pos (not long dist links) */
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL && to_pad->pos == ps->end_pos-1){
		to_pad->pos++;
		to_pad = zListMoveNext(&ps->live_to_list);
	}

	/* update and remove 3' live nodes from long dist links (pad->pos == ps->end_pos)*/
	while(to_pad != NULL && to_pad->pos == ps->end_pos){
		/* incorporate long distance live nodes (EXTERNAL state)*/
		from_pad = zListMoveFirst(&ps->live_from_list);
		while(from_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score > 
				 ps->data[from_pad->index][to_pad->state].score){				
				ps->data[from_pad->index][to_pad->state].score = 
					ps->data[from_pad->index][to_pad->index].score;
				ps->data[from_pad->index][to_pad->state].pinstate = 
					ps->data[from_pad->index][to_pad->index].pinstate;
			}
			from_pad = zListMoveNext(&ps->live_from_list);
		}
		/*remove long distance link node */
		zPinSearchResetTo(ps,to_pad->index);
		zIListAddFirst(&ps->free_to_list,to_pad->index);
		zListRemoveCurrent(&ps->live_to_list);
		to_pad = zListMoveNext(&ps->live_to_list);
	}
}

static void zStepPinSearchBack(zPinSearch* ps){
	zPinArrayData*  to_pad;
	zPinArrayData*  from_pad;

	/* Move start back */
	ps->start_pos--;

	/* create new 5' live nodes (links starting at start_pos)*/	
	zPinSearchCreateLiveNodes5(ps);

	/* foreach live to-node update path from each live from-node time start_pos */
	to_pad = zListMoveFirst(&ps->live_to_list);
	while(to_pad != NULL){
		zStepPinArrayBack(ps,to_pad);
		to_pad = zListMoveNext(&ps->live_to_list);
	}

	/* update live_from_list states pad->pos (not long dist links) */
	from_pad = zListMoveFirst(&ps->live_from_list);
	while(from_pad != NULL && from_pad->pos == ps->start_pos+1){
		from_pad->pos--;
		from_pad = zListMoveNext(&ps->live_from_list);
	}

	/* update and remove 5' live nodes from long dist links (pad->pos == ps->start_pos)*/
	while(from_pad != NULL && from_pad->pos == ps->start_pos){
		/* incorporate long distance live nodes (EXTERNAL state)*/
		to_pad = zListMoveFirst(&ps->live_to_list);
		while(to_pad != NULL){
			if(ps->data[from_pad->index][to_pad->index].score > 
				 ps->data[from_pad->state][to_pad->index].score){				
				ps->data[from_pad->state][to_pad->index].score = 
					ps->data[from_pad->index][to_pad->index].score;
				ps->data[from_pad->state][to_pad->index].pinstate = 
					ps->data[from_pad->index][to_pad->index].pinstate;
			}
			to_pad = zListMoveNext(&ps->live_to_list);
		}
		/*remove long distance link node */
		zPinSearchResetFrom(ps,from_pad->index);
		zIListAddFirst(&ps->free_from_list,from_pad->index);
		zListRemoveCurrent(&ps->live_from_list);
		from_pad = zListMoveNext(&ps->live_from_list);
	}
}

static void zStepPinArrayForward(zPinSearch* ps,zPinArrayData* from_pad){
	int            i,prestate,state;
	score_t        score;
	zIVec*         jumps;
	zHMM*          hmm = ps->trellis->hmm;

	/* reset tmp_from data */
	zPinSearchResetTmp(ps);
	
	/* Step each to state node forward */
	for(state = 0;state < ps->states;state++){
		if(hmm->state[state].type == INTERNAL ||
		   hmm->state[state].type == GINTERNAL){
			jumps = hmm->jmap[state];
			for(i = 0; i < jumps->size; i++) {
				prestate = jumps->elem[i];
				score = ps->data[from_pad->index][prestate].score;
				if(score == MIN_SCORE) continue;
				if(hmm->state[state].type == INTERNAL){
					score += zGetTransitionScore(ps->trellis->hmm,prestate,state,ps->trellis->tiso_group);
				}
				else{
					score += zGetIntergenicContinueScore(ps->trellis);
				}
				score += zScoreInternalState(ps->trellis,state,ps->end_pos,(state != prestate));
				if(score > ps->tmp[state].score){
					ps->tmp[state].score = score;
					ps->tmp[state].pinstate = ps->data[from_pad->index][prestate].pinstate;
				}
			}
		}
		else if(hmm->state[state].type == EXTERNAL){
			/* dealt with in zStepPinSearchForward */
		}
		else{
			zDie("zPinSearch only supports INTERNAL and EXTERNAL states");
		}
	}
	/* copy info from buffer back to from_pad->index */
	for(state = 0;state < ps->states;state++){
		ps->data[from_pad->index][state].score = ps->tmp[state].score;
		ps->data[from_pad->index][state].pinstate = ps->tmp[state].pinstate;
	}
}

static void zStepPinArrayBack(zPinSearch* ps,zPinArrayData* to_pad){
	int            i,poststate,state;
	score_t        score;
	zIVec*         jumps;
	zHMM*          hmm = ps->trellis->hmm;
	
	/* reset tmp_to data */
	zPinSearchResetTmp(ps);

	/* Step each from state node back */
	for(state = 0;state < ps->states;state++){
		jumps = hmm->fmap[state];
		for(i = 0; i < jumps->size; i++) {
			poststate = jumps->elem[i];
			if(hmm->state[poststate].type == INTERNAL ||
			   hmm->state[poststate].type == GINTERNAL){
				score = ps->data[poststate][to_pad->index].score;
				if(score == MIN_SCORE) continue;
				/*length 1 for duration score. it works properly since length 0 only need special
				  case for first position in sequence. see funtion in zTransition */
				if(state == 0 && poststate == 0){
					poststate = jumps->elem[i];
				}
				if(hmm->state[poststate].type == INTERNAL){
					score += zGetTransitionScore(ps->trellis->hmm,state,poststate,ps->trellis->tiso_group);
				}
				else{
					score += zGetIntergenicContinueScore(ps->trellis);
				}
				score += zScoreInternalState(ps->trellis,state,ps->start_pos,(state != poststate));
				if(score > ps->tmp[state].score){
					ps->tmp[state].score = score;
					ps->tmp[state].pinstate = ps->data[poststate][to_pad->index].pinstate;
				}
			}
			else if(hmm->state[poststate].type == EXTERNAL){
				/* dealt with when creating features (long distance links) 
					 intenrnal states preceding the external state correspinding to the 
					 feature are incorporated in zStepPinSearchBack */
			}
			else{
				zDie("zPinSearch only supports INTERNAL and EXTERNAL states");
			}
		}
	}
	/* copy info from buffer back to from_pad->index */
	for(state = 0;state < ps->states;state++){
		ps->data[state][to_pad->index].score = ps->tmp[state].score;
		ps->data[state][to_pad->index].pinstate = ps->tmp[state].pinstate;
	}
}

static void zPinSearchCreateLiveNodes5(zPinSearch* ps){
	int i;
	coor_t j;

	for (i = 0; i < ps->trellis->hmm->feature_count; i++) {
		j = ps->start_pos;
		if(ps->trellis->factory[i] != NULL){
			ps->trellis->factory[i]->create3(ps->trellis->factory[i], j, &ps->sfl);
			zPinSearchAddFeature5(ps,i,'+');
			zResetSFList(&ps->sfl);
		}
		j = ps->trellis->dna->length - ps->start_pos - 1;
		if(ps->trellis->rfactory[i] != NULL){
			ps->trellis->rfactory[i]->create5(ps->trellis->rfactory[i], j, &ps->sfl);
			zPinSearchAddFeature5(ps,i,'-');
			zResetSFList(&ps->sfl);
		}
	}	
}

static void zPinSearchCreateLiveNodes3(zPinSearch* ps){
	int i;
	coor_t j;

	for (i = 0; i < ps->trellis->hmm->feature_count; i++) {
		j = ps->end_pos;
		if(ps->trellis->factory[i] != NULL){
			ps->trellis->factory[i]->create5(ps->trellis->factory[i], j, &ps->sfl);
			zPinSearchAddFeature3(ps,i,'+');
			zResetSFList(&ps->sfl);
		}
		j = ps->trellis->dna->length - ps->end_pos - 1;
		if(ps->trellis->rfactory[i] != NULL){
			ps->trellis->rfactory[i]->create3(ps->trellis->rfactory[i], j, &ps->sfl);
			zPinSearchAddFeature3(ps,i,'-');
			zResetSFList(&ps->sfl);
		}
	}	
}

void zPinSearchAddFeature5(zPinSearch* ps, int fnum, strand_t strand){
	zTrellis*    trellis = ps->trellis;
	int          poss,i,j,k,state,prestate,poststate;
	score_t      score,base_score;
	zSfeature*   f;
	zIVec*       fmap;
	zIVec*       prejumps;
	zIVec*       postjumps;
	coor_t       tmp_start;
	zPhase_t     left_phase,right_phase;
	zPinArrayData* to_pad;
	zPinArrayData* from_pad;

	poss = ('-' != strand);
	if(poss){
		fmap = trellis->fmap3[fnum];
	}
	else{
		fmap = trellis->rfmap3[fnum];
	}
	
	f = zSFListMoveFirst(&ps->sfl);
	/* add new live node for each possible start state of each feature in sfl */
	while(f != NULL){
		f->strand = strand;	
		if(poss){
			right_phase = zSFEndPhase(trellis->dna, f);
			left_phase = zSFBeginPhase(trellis->dna, f);
		}
		else{
			tmp_start = f->start;
			f->start  = trellis->dna->length - f->end - 1;
			f->end    = trellis->dna->length - tmp_start - 1;			
			right_phase = zSFBeginPhase(trellis->rdna, f);
			left_phase = zSFEndPhase(trellis->rdna, f);
		}
		if(f->end != ps->start_pos){
			zDie("Cannot add external features to 5' end of pin search which don't end at start_pos");
		}
		/* for each (external) state which this feature can end in */
		for(i = 0;i < fmap->size;i++){
			state = fmap->elem[i];
			/* make sure external state (state) is compatible with feature */
			if(right_phase != trellis->hmm->state[state].phase) continue;
			/* for each (internal) state which can transition into this (external) state */ 
			prejumps = trellis->hmm->jmap[state];
			for(j = 0; j < prejumps->size;j++){
				prestate = prejumps->elem[j];
				if (!zCompatPhases(poss, trellis->hmm->state[prestate].phase, left_phase)) 
					continue;
				/* add a new from pad */
				from_pad = zPinSearchGetFromPad(ps,f->start-1,prestate);
				from_pad->link_state = state;
				/* for each (internal) state which this (external) state can transition to */
				postjumps = trellis->hmm->fmap[state];
				for(k = 0; k < postjumps->size;k++){
					poststate = postjumps->elem[k];
					base_score = zScoreInternalState(trellis,poststate,f->end+1,true) +
						zGetTransitionScore(ps->trellis->hmm,state,poststate,ps->trellis->tiso_group)+
						zScoreExternalState(trellis,state,f->end,f->end-f->start+1,f->score) + 
						zGetTransitionScore(ps->trellis->hmm,prestate,state,ps->trellis->tiso_group);
					
					/* select the best route to each to_pad.  each route can take any
						 (internal) prestate -> (external) state -> (internal) poststate */
					/*fprintf(stderr,"LDL %d -> %d (%d) -> %d\n",prestate,state,fnum,poststate);EVAN*/
					to_pad = zListMoveFirst(&ps->live_to_list);
					while(to_pad != NULL){
						score = base_score + ps->data[poststate][to_pad->index].score;
						if(score > ps->data[from_pad->index][to_pad->index].score){
							ps->data[from_pad->index][to_pad->index].score = score;
							ps->data[from_pad->index][to_pad->index].pinstate = 
								ps->data[poststate][to_pad->index].pinstate;
						}
						to_pad = zListMoveNext(&ps->live_to_list);
					}
				}
			}
		}
		f = zSFListMoveNext(&ps->sfl);
	}
}

void zPinSearchAddFeature3(zPinSearch* ps, int fnum, strand_t strand){
	zTrellis*      trellis = ps->trellis;
	int            poss,i,j,state,prestate;
	score_t        score,base_score;
	zSfeature*     f;
	zIVec*         fmap;
	zIVec*         jumps;
	coor_t         tmp_start;
	zPhase_t       left_phase,right_phase;
	zPinArrayData* to_pad;
	zPinArrayData* from_pad;
	
	poss = ('-' != strand);
	if(poss){
		fmap = trellis->fmap3[fnum];
	}
	else{
		fmap = trellis->rfmap3[fnum];
	}

	f = zSFListMoveFirst(&ps->sfl);
	/* add new live node for each state reachable by each feature in sfl */
	while(f != NULL){
		f->strand = strand;	
		if(poss){
			right_phase = zSFEndPhase(trellis->dna, f);
			left_phase = zSFBeginPhase(trellis->dna, f);
		}
		else{
			tmp_start = f->start;
			f->start  = trellis->dna->length - f->end - 1;
			f->end    = trellis->dna->length - tmp_start - 1;			
			right_phase = zSFBeginPhase(trellis->rdna, f);
			left_phase = zSFEndPhase(trellis->rdna, f);
		}			
		if(f->start != ps->end_pos){
			zDie("Cannot add external features to 3' end of pin search which don't begin at end_pos");
		}
		/* for each (external) state which this feature can end in */
		for(i = 0;i < fmap->size;i++){
			state = fmap->elem[i];
			/* make sure extenral state (state) is compatible with feature */
			if(right_phase != trellis->hmm->state[state].phase) continue;					
			to_pad = zPinSearchGetToPad(ps,f->end,state);
			/* for each from_state select the best path to f->end in state to_pad->state.  
				 This path may go through any of the states which transition to to_pad->state */
			jumps = trellis->hmm->jmap[state];
			for(j = 0;j < jumps->size;j++){
				prestate = jumps->elem[j];
				/* make sure internal state (prestate) is compatible with feature */
				if (!zCompatPhases(poss, trellis->hmm->state[prestate].phase, left_phase)) 
					continue;
				/* for each possible state state/pos, if path to end state though prestate 
					 is the best we've seen so far use it */
				base_score = zScoreExternalState(trellis,state,f->end,f->end-f->start+1,
												 f->score) +
					zGetTransitionScore(trellis->hmm,prestate,state,trellis->tiso_group);
				from_pad = zListMoveFirst(&ps->live_from_list);
				while(from_pad != NULL){
					score = base_score + ps->data[from_pad->index][prestate].score;
					if(score > ps->data[from_pad->index][to_pad->index].score){
						ps->data[from_pad->index][to_pad->index].score = score;
						ps->data[from_pad->index][to_pad->index].pinstate = 
							ps->data[from_pad->index][prestate].pinstate;
					}
					from_pad = zListMoveNext(&ps->live_from_list);
				}
			}
		}
		f = zSFListMoveNext(&ps->sfl);
	}
}

void zPinSearchAddFeatureInit(zPinSearch* ps, int fnum, strand_t strand){
	zTrellis*      trellis = ps->trellis;
	int            poss,i,j,state,prestate;
	score_t        score;
	zSfeature*     f;
	zIVec*         fmap;
	zIVec*         jumps;
	coor_t         tmp_start;
	zPhase_t       left_phase,right_phase;
	zPinArrayData* to_pad;
	zPinArrayData* from_pad;
	
	poss = ('-' != strand);
	if(poss){
		fmap = trellis->fmap3[fnum];
	}
	else{
		fmap = trellis->rfmap3[fnum];
	}

	f = zSFListMoveFirst(&ps->sfl);
	/* add new live to and from nodes for each state reachable by each feature in sfl */
	while(f != NULL){

		zCheckPinSearchListOrder(ps);

		/*fprintf(stderr,"Add %d\t%d -> %d %c\n",fnum,f->start,f->end,f->strand);EVAN*/
		f->strand = strand;	
		if(poss){
			right_phase = zSFEndPhase(trellis->dna, f);
			left_phase = zSFBeginPhase(trellis->dna, f);
		}
		else{
			tmp_start = f->start;
			f->start  = trellis->dna->length - f->end - 1;
			f->end    = trellis->dna->length - tmp_start - 1;			
			right_phase = zSFBeginPhase(trellis->rdna, f);
			left_phase = zSFEndPhase(trellis->rdna, f);
		}			
		/* we need at least one pos before f->start for the (internal) prestate
			 this will never be a feature that would be used anyway and no search 
			 will be successful if started this close to the end of the sequence */
		if(f->start == 0) continue;
		/* for each (external) state which this feature can end in */
		for(i = 0;i < fmap->size;i++){
			state = fmap->elem[i];
			/* make sure extenral state (state) is compatible with feature */
			if(right_phase != trellis->hmm->state[state].phase) continue;					
			to_pad = zPinSearchGetToPad(ps,f->end,state);
			
			/* for each state which transitions to to_pad->state */
			jumps = trellis->hmm->jmap[state];
			for(j = 0;j < jumps->size;j++){
				prestate = jumps->elem[j];
				/* make sure internal state (prestate) is compatible with feature */
				if (!zCompatPhases(poss, trellis->hmm->state[prestate].phase, left_phase)) 
					continue;
				/* for each possible state state/pos, if path to end state though prestate 
				   is the best we've seen so far use it */
				score = zScoreInternalState(trellis,prestate,f->start-1,true)+
					zGetTransitionScore(trellis->hmm,prestate,state,trellis->tiso_group) +
					zScoreExternalState(trellis,state,f->end,f->end-f->start+1,
										f->score);
				from_pad = zPinSearchGetFromPad(ps,f->start-1,prestate);
				if(score > ps->data[from_pad->index][to_pad->index].score){
					ps->data[from_pad->index][to_pad->index].score = score;
					ps->data[from_pad->index][to_pad->index].pinstate = state;
				}			
			}
		}
		f = zSFListMoveNext(&ps->sfl);
	}
}

#endif

