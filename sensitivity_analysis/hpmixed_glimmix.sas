

%let dir = C:\DATASETS; 

PROC IMPORT OUT=dat DATAFILE="&dir\dat_gee.csv" REPLACE dbms=csv; GETNAMES=yes; RUN; 

%macro enrichment (annot, sig);
proc hpmixed data = dat; class seqnames; 
model &annot = &sig; 
random seqnames; 
ods output CovParms=HPMEstimate; run;

proc glimmix data=dat; 
class seqnames; 
model &annot = &sig / dist = binomial link = logit s; 
random int /subject = seqnames; 
parms /pdata=HPMEstimate hold=1,2 noiter;
ods output ParameterEstimates = result; 
run;

%let effect = "&sig";
data pvalsOR; set result; 
if Effect = &effect; 
or = exp(Estimate); run;

proc print data = pvalsOR;
format probt e10.;
run;

proc freq data = dat; 
tables &sig*&annot/ nocol nopct or; run;
%mend;

*supp3 -- hyper;
%enrichment(isIsland,isSig_hyper);
%enrichment(isOpenSea,isSig_hyper);
%enrichment(isShelf,isSig_hyper);
%enrichment(isShore,isSig_hyper);
%enrichment(is1stExon,isSig_hyper);
%enrichment(is3UTR,isSig_hyper);
%enrichment(is5UTR,isSig_hyper);
%enrichment(isBody,isSig_hyper);
%enrichment(isIntergenic,isSig_hyper);
%enrichment(isTSS1500,isSig_hyper);
%enrichment(isTSS200,isSig_hyper);

*supp3 -- hypo;
%enrichment(isIsland,isSig_hypo);
%enrichment(isOpenSea,isSig_hypo);
%enrichment(isShelf,isSig_hypo);
%enrichment(isShore,isSig_hypo);
%enrichment(is1stExon,isSig_hypo);
%enrichment(is3UTR,isSig_hypo);
%enrichment(is5UTR,isSig_hypo);
%enrichment(isBody,isSig_hypo);
%enrichment(isIntergenic,isSig_hypo);
%enrichment(isTSS1500,isSig_hypo);
%enrichment(isTSS200,isSig_hypo);

*supp4 -- hyper;
%enrichment(isActiveTSS,isSig_hyper);
%enrichment(isBivalentEnhancer,isSig_hyper);
%enrichment(isBivalentPoisedTSS,isSig_hyper);
%enrichment(isEnhancers,isSig_hyper);
%enrichment(isFlankingActiveTSS,isSig_hyper);
%enrichment(isFlankingBivalentTSSEnh,isSig_hyper);
%enrichment(isGenicenhancers,isSig_hyper);
%enrichment(isHeterochromatin,isSig_hyper);
%enrichment(isQuiescentLow,isSig_hyper);
%enrichment(isRepressedPolyComb,isSig_hyper);
%enrichment(isStrongtranscription,isSig_hyper);
%enrichment(isTranscratgene5and3,isSig_hyper);
%enrichment(isWeakRepressedPolyComb,isSig_hyper);
%enrichment(isWeaktranscription,isSig_hyper);
%enrichment(isZNFgenesrepeats,isSig_hyper);

*supp4 -- hypo;
%enrichment(isActiveTSS,isSig_hypo);
%enrichment(isBivalentEnhancer,isSig_hypo);
%enrichment(isBivalentPoisedTSS,isSig_hypo);
%enrichment(isEnhancers,isSig_hypo);
%enrichment(isFlankingActiveTSS,isSig_hypo);
%enrichment(isFlankingBivalentTSSEnh,isSig_hypo);
%enrichment(isGenicenhancers,isSig_hypo);
%enrichment(isHeterochromatin,isSig_hypo);
%enrichment(isQuiescentLow,isSig_hypo);
%enrichment(isRepressedPolyComb,isSig_hypo);
%enrichment(isStrongtranscription,isSig_hypo);
%enrichment(isTranscratgene5and3,isSig_hypo);
%enrichment(isWeakRepressedPolyComb,isSig_hypo);
%enrichment(isWeaktranscription,isSig_hypo);
%enrichment(isZNFgenesrepeats,isSig_hypo);

