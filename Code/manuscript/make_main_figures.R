# =============================================================================
# make_main_figures.R -- main-text Figures 1-4 of the manuscript and the sidecar
# CSVs that record every plotted point and interval. Reads accepted saved results
# only: no topic-model fitting, no rescoring. Every input is checked (internal
# stopifnot gates) against the accepted revision outputs before plotting.
#   Figure 1  F1_selection_support.pdf   (fig:e1_curves)
#   Figure 2  gap_decomposition.pdf      (fig:gap_channels)
#   Figure 3  F2_mdna_fit.pdf            (fig:mdna_fit)
#   Figure 4  residual_mass.pdf          (fig:mdna_mass)
# Inputs : Data/E1/e1_results_..._rev2.qs2, Data/MDNA/mdna_results_MDNA_2015_2016_rev2.qs2,
#          Data/MDNA/mdna_prep_2015_2016.qs2 (Zenodo results record), Results/csv/*.
# Outputs: Results/manuscript/ (4 PDFs + 13 CSV sidecars).
# Usage  : Rscript Code/manuscript/make_main_figures.R      (from the package root; ~5 s)
# Origin : the authors' working copy make_short_figures.R (24 Sep 2026); only the
#          root and output paths differ.
# =============================================================================
# Plot and summarise accepted saved scores only: no topic-model fitting or rescoring.
suppressPackageStartupMessages({library(qs2);library(data.table);library(ggplot2);library(patchwork)})
root <- here::here()
e1 <- qs_read(file.path(root,'Data/E1/e1_results_E1_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01_rev2.qs2'))
m <- qs_read(file.path(root,'Data/MDNA/mdna_results_MDNA_2015_2016_rev2.qs2'))
u <- fread(file.path(root,'Results/csv/rev2_unbinned_comparator_rev2.csv'))
g <- fread(file.path(root,'Results/csv/mdna_gain_profile_rev1.csv'))[protocol=='com']
outdir <- file.path(root,'Results/manuscript'); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
labs <- c(insample='In-sample',ho_reconstruction='Reconstruction',ho_completion='Completion',ins='In-sample',rec='Reconstruction',com='Completion')
protocol_factor <- function(x) factor(labs[x],levels=c('In-sample','Reconstruction','Completion'))
colours <- c('In-sample'='#777777','Reconstruction'='#0072B2','Completion'='#B34C00')
line_types <- c('In-sample'='solid','Reconstruction'='dashed','Completion'='longdash')
protocol_scales <- function() list(
  scale_colour_manual(values=colours,limits=names(colours)),
  scale_fill_manual(values=colours,limits=names(colours)),
  scale_linetype_manual(values=line_types,limits=names(colours)),
  guides(fill='none',linetype='none',colour=guide_legend(override.aes=list(linetype=unname(line_types)))))
style <- theme_bw(base_size=11,base_family='serif')+
  theme(panel.grid.minor=element_blank(),legend.position='bottom',legend.title=element_blank(),
        legend.text=element_text(size=9),legend.key.width=grid::unit(.9,'cm'),
        plot.title=element_text(size=11),plot.margin=margin(5,9,5,5))
# Figure 1A: the unit of Monte Carlo replication is one independent training
# seed with its evaluation corpus, not one document or one nested evaluation.
seed <- e1$summary[metric=='dev',.(K,eval,replicate,r2_macro)]
stopifnot(!anyDuplicated(seed[,.(K,eval,replicate)]))
s <- seed[,.(fit=mean(r2_macro),seed_sd=sd(r2_macro),n_seed=.N),by=.(K,eval)]
stopifnot(nrow(s)==30L,all(s$n_seed==10L))
s[,`:=`(mc_se=seed_sd/sqrt(n_seed),t_critical=qt(.975,n_seed-1L))]
s[,`:=`(lwr=fit-t_critical*mc_se,upr=fit+t_critical*mc_se,Protocol=protocol_factor(eval))]
p1 <- ggplot(s,aes(K,fit,linetype=Protocol,colour=Protocol))+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Protocol),alpha=.14,colour=NA,show.legend=FALSE)+
  geom_vline(xintercept=40,colour='grey65')+geom_line(linewidth=.65)+
  protocol_scales()+
  labs(x='Number of topics, K',y='Deviance Macro',title='A. Ten-run mean and 95% MC intervals')+style
b <- e1$summary[metric=='dev' & eval=='ho_completion' & replicate==1L,.(K,fit=r2_micro,Support='Harmonised')]
a <- u[source=='E1_Kstar40_W5000' & protocol=='com' & floor==1e-12,.(K,fit=r2_micro,Support='Unbinned')]
stopifnot(nrow(a)==10L,nrow(b)==10L)
v <- rbind(b,a)
p2 <- ggplot(v,aes(K,fit,linetype=Support))+geom_vline(xintercept=40,colour='grey65')+
  geom_line(linewidth=.65)+labs(x='Number of topics, K',y='Completion Deviance Micro',title='B. Support comparison: seed 1')+style
# Figure 2A: pointwise held-out intervals conditional on the empirical training
# fit, under working document independence. No inferential band for in-sample fit.
t <- copy(m$summary[metric=='dev'])
t[,Protocol:=protocol_factor(eval)]
heldout <- rbindlist(list(copy(m$doc_rec)[,eval:='rec'],copy(m$doc_com)[,eval:='com']))[
  metric=='dev' & is.finite(r2_doc) & d_null>0 & d_null>=1]
ci <- heldout[,.(mean_check=mean(r2_doc),n=.N,sd=sd(r2_doc),se=sd(r2_doc)/sqrt(.N)),by=.(K,eval)]
t <- merge(t,ci,by=c('K','eval'),all.x=TRUE,sort=FALSE)
stopifnot(max(abs(t[eval!='ins',r2_macro-mean_check]))<1e-10)
t[,`:=`(lwr=r2_macro-qnorm(.975)*se,upr=r2_macro+qnorm(.975)*se)]
old_ci <- m$battery_ci[metric=='dev' & eval%in%c('rec','com'),.(K,eval,lwr_saved=lwr,upr_saved=upr)]
ci_check <- merge(t,old_ci,by=c('K','eval'))
stopifnot(nrow(ci_check)>0L,max(abs(ci_check$lwr-ci_check$lwr_saved))<1e-10,
          max(abs(ci_check$upr-ci_check$upr_saved))<1e-10)

# Independent transcription check of the entire 190-pair family, calculated
# from stored document scores and compared with the accepted gain-profile CSV.
w <- dcast(heldout[eval=='com'],doc_id~K,value.var='r2_doc')
ks <- sort(unique(heldout[eval=='com',K]));stopifnot(length(ks)==20L,nrow(w)==3447L,!anyNA(w))
M <- choose(length(ks),2); pairs <- rbindlist(lapply(seq_len(length(ks)-1),function(i)
  rbindlist(lapply((i+1):length(ks),function(j){
    d <- w[[as.character(ks[j])]]-w[[as.character(ks[i])]]
    data.table(K=ks[i],K_to=ks[j],gain=mean(d),se=sd(d)/sqrt(length(d)))
  }))))
pairs[,ub:=gain+qnorm(1-.05/M)*se]
recalc <- pairs[,.(mean_gain_check=max(gain),upper_check=max(ub)),by=K]
g <- merge(g,recalc,by='K')
stopifnot(nrow(pairs)==190L,nrow(g)==19L,
          max(abs(g$mean_gain_check-g$total_gain_max))<1e-10,
          max(abs(g$upper_check-g$total_gain_ub))<1e-10,
          g[total_gain_ub<=.01,min(K)]==180L,g[total_gain_ub<=.005,min(K)]==190L)
pa <- ggplot(t,aes(K,r2_macro,linetype=Protocol,colour=Protocol))+
  geom_ribbon(data=t[eval!='ins'],aes(ymin=lwr,ymax=upr,fill=Protocol),alpha=.16,colour=NA,show.legend=FALSE)+
  geom_line(linewidth=.65)+protocol_scales()+
  labs(x='Number of topics, K',y='Deviance Macro',title='A. Fit and pointwise 95% intervals')+style
pb <- ggplot(g,aes(K,total_gain_ub))+geom_line(linewidth=.65)+geom_point(size=1)+
  geom_hline(yintercept=.01,linetype='dashed',colour='#0072B2')+
  geom_hline(yintercept=.005,linetype='dotted',colour='#B34C00')+
  annotate('text',x=15,y=.0115,label='Tolerance: 0.01',hjust=0,vjust=0,size=3,family='serif',colour='#0072B2')+
  annotate('text',x=15,y=.0057,label='Tolerance: 0.005',hjust=0,vjust=0,size=3,family='serif',colour='#B34C00')+
  geom_point(data=g[K%in%c(180L,190L)],size=2)+
  annotate('text',x=163,y=.0074,label='180',size=3,family='serif')+
  annotate('text',x=170,y=.0024,label='190',size=3,family='serif')+
  scale_y_log10(breaks=c(.002,.005,.01,.03,.1,.3),labels=c('0.002','0.005','0.01','0.03','0.1','0.3'))+
  scale_x_continuous(limits=c(10,200),breaks=c(20,60,100,140,180))+
  labs(x='Number of topics, K',y='Upper bound on remaining gain',title='B. Completion: simultaneous total gain')+style
# Save sidecars for every plotted point and interval, including the 10 raw runs.
fwrite(seed,file.path(outdir,'F1_seed_scores.csv'))
fwrite(s,file.path(outdir,'F1_panel_A.csv'));fwrite(v,file.path(outdir,'F1_panel_B.csv'))
fwrite(t,file.path(outdir,'F2_data.csv'));fwrite(g,file.path(outdir,'F2_panel_B.csv'))
fwrite(pairs,file.path(outdir,'F2_gain_pairs.csv'))
ggsave(file.path(outdir,'F1_selection_support.pdf'),p1+p2,width=6.8,height=3.0,device=cairo_pdf)
ggsave(file.path(outdir,'F2_mdna_fit.pdf'),pa+pb,width=6.8,height=3.05,device=cairo_pdf)
cat('Figure checks passed: 30 ten-run intervals; 40 held-out empirical intervals; 190 paired gains; selections 180 and 190.\n')

# Aggregation identity: mean components, with MC uncertainty for the total gap.
e3 <- fread(file.path(root,'Results/csv/e3_gap_decomposition_E3_full_Kstar40_J1000_W5000_a0p5_b0p01_k10-100by10_warplda_fa0p1_fb0p01.csv'))[
  K==40 & eval=='ho_reconstruction']
stopifnot(nrow(e3)==30L,!anyDuplicated(e3[,.(scenario,replicate)]),
          max(abs(e3$gap-e3$ch_length-e3$ch_atypicality-e3$ch_interaction))<1e-12)
eg <- e3[,.(gap=mean(gap),ch_length=mean(ch_length),ch_atypicality=mean(ch_atypicality),
             ch_interaction=mean(ch_interaction),n_seed=.N,mc_se=sd(gap)/sqrt(.N)),by=scenario]
eg[,`:=`(lwr=gap-qt(.975,9)*mc_se,upr=gap+qt(.975,9)*mc_se)]
scenario_labels <- c(A='A: heterogeneous lengths',B='B: approximately equal lengths',C='C: fixed length, mixed concentration')
eg[,Scenario:=factor(scenario_labels[scenario],levels=rev(scenario_labels))]
egl <- melt(eg,id.vars=c('scenario','Scenario'),measure.vars=c('ch_length','ch_atypicality','ch_interaction'),variable.name='component',value.name='contribution')
egl[,Component:=factor(component,levels=c('ch_length','ch_atypicality','ch_interaction'),labels=c('Length','Atypicality','Interaction'))]
pgap <- ggplot(egl,aes(contribution,Scenario,fill=Component))+
  geom_col(width=.55,position=position_stack(reverse=TRUE),colour='white',linewidth=.15)+
  geom_errorbar(data=eg,aes(xmin=lwr,xmax=upr,y=Scenario),inherit.aes=FALSE,width=.18,orientation='y',linewidth=.45)+
  geom_point(data=eg,aes(x=gap,y=Scenario),inherit.aes=FALSE,shape=18,size=2.6)+
  scale_fill_manual(values=c(Length='#d4d4d4',Atypicality='#929292',Interaction='#4d4d4d'))+
  scale_x_continuous(breaks=seq(0,.2,.05),limits=c(0,.215),expand=expansion(mult=c(0,.02)))+
  labs(x='Deviance Micro - Macro',y=NULL)+style+
  theme(legend.key.width=grid::unit(.45,'cm'),panel.grid.major.y=element_blank())
fwrite(eg,file.path(outdir,'gap_decomposition_data.csv'))
fwrite(e3,file.path(outdir,'gap_decomposition_seed_data.csv'))
ggsave(file.path(outdir,'gap_decomposition.pdf'),pgap,width=6.8,height=2.6,device=cairo_pdf)

# Figure 4: fixed training-frequency groups, pairing per-word contrasts and
# group totals. Both panels use firm-clustered pointwise 95% normal intervals.
rmass <- fread(file.path(root,'Results/csv/rev2_mdna_residual_mass_rev2.csv'))[
  K%in%c(50L,180L) & partition=='training frequency (Test 2)']
contrast <- fread(file.path(root,'Results/csv/rev2_mdna_moment_strata_rev2.csv'))[
  K%in%c(50L,180L) & test=='T2_freq_strata']
stopifnot(nrow(rmass)==10L,nrow(contrast)==8L,all(is.finite(contrast$t)),all(contrast$t!=0))
contrast[, `:=`(se=abs(gbar/t), group=as.integer(sub('f([1-4])_vs_f5','\\1',stratum)))]
contrast[, se_iid:=se]
setorder(contrast,K,group); setorder(rmass,K,stratum)
# Independently link the two panels: difference of mass per word equals contrast.
for (i in seq_len(nrow(contrast))) {
  a <- contrast[i]; lo <- rmass[K==a$K & stratum==a$group]; hi <- rmass[K==a$K & stratum==5L]
  stopifnot(nrow(lo)==1L,nrow(hi)==1L,
            abs(a$gbar-(lo$mass_pp/(100*lo$words)-hi$mass_pp/(100*hi$words)))<1e-15)
}
# Recover document-level group masses from the four cached contrasts.
# If u_b is mean residual per word and n_b group size, g_b=u_b-u_5;
# sum_b n_b*u_b=0 implies u_5=-sum_{b<5} n_b*g_b/W.
# This uses the exhaustive frequency partition, not an assumption of zero bias.
prep_path <- file.path(root,'Data/MDNA/mdna_prep_2015_2016.qs2')
prep <- qs_read(prep_path)
stopifnot(identical(digest::digest(prep$dtm_ev,algo='xxhash64'),m$input_check$dtm_ev_hash),
          identical(digest::digest(prep$dtm_train,algo='xxhash64'),m$input_check$dtm_train_hash))
source(file.path(root,'Code/R/utils_moment_tests.R'))
ns <- tabulate(.freq_strata(prep$dtm_train,B=5L),nbins=5L)
stopifnot(sum(ns)==ncol(prep$dtm_train),all(ns==rmass[K==50L][order(stratum),words]),
          all(ns==rmass[K==180L][order(stratum),words]))
cluster_rows <- list(); meta_rows <- list()
crit <- qnorm(.975)
accepted_tests <- fread(file.path(root,'Results/csv/rev2_mdna_moment_tests_rev2.csv'))
for (kv in c(50L,180L)) {
  G <- m$moments_ho[[as.character(kv)]]$T2_freq_strata
  stopifnot(identical(rownames(G),rownames(prep$dtm_ev)),
            identical(colnames(G),paste0('f',1:4,'_vs_f5')))
  idx <- match(rownames(G),prep$dv_ev$doc_id)
  cl <- as.character(prep$dv_ev$cik[idx]); stopifnot(!anyNA(cl),!anyDuplicated(rownames(G)))
  ref <- -as.vector(G %*% ns[1:4])/sum(ns)
  R <- sweep(cbind(sweep(G,1L,ref,'+'),ref),2L,ns,'*')*100
  stopifnot(max(abs(rowSums(R)))<1e-12,
            max(abs(colMeans(R)-rmass[K==kv][order(stratum),mass_pp]))<1e-10,
            max(abs(colMeans(G)-contrast[K==kv][order(group),gbar]))<1e-15)
  Y <- cbind(G,R); colnames(Y) <- c(paste0('c',1:4),paste0('m',1:5))
  centred <- sweep(Y,2L,colMeans(Y))
  S <- rowsum(centred,cl); J <- nrow(G); C <- nrow(S)
  stopifnot(J==3466L,C==3059L)
  V <- C/(C-1)*crossprod(S)/J^2
  cluster_se <- sqrt(diag(V)); mu <- colMeans(G)
  stat <- as.numeric(t(mu)%*%solve(V[1:4,1:4])%*%mu)
  stopifnot(abs(stat-accepted_tests[K==kv & test=='T2_freq_strata',stat_cluster])<1e-8)
  contrast[K==kv,`:=`(se=cluster_se[1:4],lwr=gbar-crit*cluster_se[1:4],upr=gbar+crit*cluster_se[1:4],
                     pval_cluster=2*pnorm(-abs(gbar/cluster_se[1:4])),n_docs=J,n_clusters=C,critical=crit)]
  rmass[K==kv,`:=`(se=cluster_se[5:9],lwr=mass_pp-crit*cluster_se[5:9],upr=mass_pp+crit*cluster_se[5:9],
                  n_docs=J,n_clusters=C,critical=crit)]
  # Anonymous centred cluster sums suffice to independently audit every interval.
  cluster_rows[[as.character(kv)]] <- cbind(data.table(K=kv,cluster_index=seq_len(C)),as.data.table(S))
  meta_rows[[as.character(kv)]] <- data.table(K=kv,n_docs=J,n_clusters=C,critical=crit,
                                            correction=C/(C-1),joint_stat=stat)
}
fwrite(rbindlist(cluster_rows),file.path(outdir,'residual_cluster_sums.csv'))
fwrite(rbindlist(meta_rows),file.path(outdir,'residual_interval_metadata.csv'))
fwrite(data.table(source=c('Data/MDNA/mdna_results_MDNA_2015_2016_rev2.qs2',
                          'Data/MDNA/mdna_prep_2015_2016.qs2'),
                 md5=unname(tools::md5sum(c(file.path(root,'Data/MDNA/mdna_results_MDNA_2015_2016_rev2.qs2'),prep_path)))),
       file.path(outdir,'residual_interval_sources.csv'))
chk <- rmass[,.(mass_sum=sum(mass_pp),half_check=sum(abs(mass_pp))/2,
                half_saved=unique(half_abs_sum_pp)),by=K]
stopifnot(max(abs(chk$mass_sum))<1e-9,max(abs(chk$half_check-chk$half_saved))<1e-12)
for (d in list(rmass,contrast)) d[,Model:=factor(paste('K =',K),levels=c('K = 50','K = 180'))]
shape_scale <- function() scale_shape_manual(values=c(16,2))
pcontrast <- ggplot(contrast,aes(factor(group),gbar*1e7,shape=Model,group=Model))+
  geom_hline(yintercept=0,colour='grey45',linewidth=.4)+
  geom_errorbar(aes(ymin=lwr*1e7,ymax=upr*1e7),position=position_dodge(width=.32),
                width=.16,linewidth=.45)+
  geom_point(position=position_dodge(width=.32),size=2.5,stroke=.7)+shape_scale()+
  scale_x_discrete(labels=paste(1:4,'vs 5'))+
  labs(title='A. Direction and uncertainty',x='Frequency group versus highest',
       y=expression(atop('Mean residual contrast',paste('(',10^{-7},' probability per word)'))))+style
prmass <- ggplot(rmass,aes(factor(stratum),mass_pp,shape=Model,group=Model))+
  geom_hline(yintercept=0,colour='grey45',linewidth=.4)+
  geom_errorbar(aes(ymin=lwr,ymax=upr),position=position_dodge(width=.32),width=.16,linewidth=.45)+
  geom_point(position=position_dodge(width=.32),size=2.5,stroke=.7)+shape_scale()+
  scale_y_continuous(breaks=c(-.08,-.04,0,.04,.08))+
  labs(title='B. Aggregate magnitude',x='Frequency group (lowest to highest)',
       y='Residual mass\n(percentage points)')+style
fig4 <- (pcontrast+prmass+plot_layout(guides='collect')) &
  theme(legend.position='bottom',legend.key.width=grid::unit(.5,'cm'))
fwrite(rmass,file.path(outdir,'residual_mass_data.csv'))
fwrite(contrast,file.path(outdir,'residual_contrast_data.csv'))
ggsave(file.path(outdir,'residual_mass.pdf'),fig4,width=6.8,height=2.8,device=cairo_pdf)
cat('Figure checks passed: 30 seed decompositions; 18 firm-clustered intervals; unchanged frequency contrasts and group masses; cross-panel identity.\n')
