% Demo with AFNI data

M=spm_read_vols(spm_vol('full_mask.FT.nii.gz'));
T=spm_read_vols(spm_vol('stats.FT.nii.gz'));
T=T(M>0);
df=430;
alph=0.05;
P=tcdf(T,df,'upper');
P2=2*tcdf(abs(T),df,'upper');
Ipos=T>0;
Ineg=T<0;
Iref=@(x)(1:length(x))/length(x);
xx=-10:0.01:10;

figure
subplot(2,1,1)
qt=fdr(P2);
tt=tinv(1-qt/2,df);
histogram(T,100,'Normalization','pdf');xline(-tt);xline(tt);
hold on; line(xx,tpdf(xx,df),'color','green','linewidth',2); hold off
title({sprintf('Two-sided FDR: Thr P=%f T=%.3f,%.3f',qt,-tt,tt),
       sprintf('%d voxels detected',sum(P2<=qt))})
subplot(2,1,2)
loglog(Iref(P2),sort(P2),'.');grid on;refline(1,0);refline(alph,0);title('Two-sided FDR')

figure
q1=fdr(P);
t1=tinv(1-q1,df);
q2=fdr(1-P);
t2=-tinv(1-q2,df);
subplot(2,1,1)
histogram(T,100,'Normalization','pdf');xline(t1);xline(t2)
hold on; line(xx,tpdf(xx,df),'color','green','linewidth',2); hold off
title({sprintf('FDR Twice: Thr P=%f,%f T=%.3f,%.3f',q2,q1,t2,t1),
       sprintf('%d voxels detected',sum(P<=q1)+sum(1-P<=q2))})
subplot(2,1,2)
loglog(Iref(P),sort(P),'.');
hold on
loglog(Iref(1-P),sort(1-P),'.');
hold off
grid on;abline(0,1,'linestyle','-');abline(0,alph);title('FDR Twice')

figure
q1p=fdr(P2(Ipos));
q2p=fdr(P2(Ineg));
t1p=tinv(q1p/2,df);
t2p=-tinv(q2p/2,df);
subplot(2,1,1)
histogram(T,100,'Normalization','pdf');xline(t1p);xline(t2p)
hold on; line(xx,tpdf(xx,df),'color','green','linewidth',2); hold off
title({sprintf('FDR split pos/neg: Thr P=%f,%f T=%.3f,%.3f',q2p,q1p,t2p,t1p),
       sprintf('%d voxels detected',sum(P<=q1p)+sum(1-P<=q2p))})
subplot(2,1,2)
loglog(Iref(P2(Ipos)),sort(P2(Ipos)),'.');
hold on
loglog(Iref(P2(Ineg)),sort(P2(Ineg)),'.');
hold off
grid on;abline(0,1,'linestyle','-');abline(0,alph);title('FDR split pos/neg')

