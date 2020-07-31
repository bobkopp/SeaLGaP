function [gsl,iced,steric] = SampleLR04BasedPrior(Nsamp,icedistfile,datset,t2,ctrtime,ctrshiftsd,timestepsd,scal,scalsd,timescale,esl0)

% [gsl,iced,steric] = SampleLR04BasedPrior(Nsamp,icedistfile,datset,t2,ctrtime,ctrshiftsd,timestepsd,scal,scalsd,timescale,esl0)
%
% Sample prior based on Lisiecki-Raymo global oxygen isotope curve.
%
% Last updated by Bob Kopp rkopp-at-princeton.edu, 13 July 2009

defval('datset','LR04_LIG');
defval('t2',150:-1:90);
defval('Nsamp',500);

defval('icedistfile','icedist');

defval('ctrtime',120);
defval('ctrshiftsd',0);
defval('timestepsd',.25);

% assume LGM sea level of 125 m with change in d18O of 1.79 per mil yields a slope of 70 m/permil
% from Bintanja et al., we assume a 2 sigma uncertainty of +/- 20 m in the projection
% (sd of 10 m)

defval('scal',70);
defval('projivsd',10);

% and we'll put a fairly weak constraint on the slope as used for the derivative
% from Bintanja et al, it looks like the slope can range from 25% to 150% of 
% mean value, so we'll use a 2 sigma slope uncertainty of +/- 50 m/permil
% (sd of 25 m/permil)

defval('scalsd',25)

defval('timescale',3); % this is the timescale over which scal is changing

defval('esl0',[0.18 7.3 0.02 5.0 58.0 0 0]);

defval('stericsd',2);
defval('sterictimescale',2);

LR04_Present = [0.0 3.23 0.03];


if isstr(icedistfile)
	load(icedistfile);
elseif iscell(icedistfile)
	mu=icedistfile{1};
	cv=icedistfile{2};
	bins=icedistfile{3};
end

LR04=LookupLR04(datset);

A=meshgrid(LR04(:,2));
Asd=meshgrid(LR04(:,3));
d18Obase=LR04(:,2)-LR04_Present(2);
d18Obasesd=sqrt(LR04(:,3).^2+LR04_Present(3)^2); d18Obasesd(1) = LR04_Present(3);


ivbase=d18Obase*scal;
ivsd=ones(size(ivbase))*projivsd;

diffd18O=A-A';
diffd18Osd=sqrt(Asd.^2+Asd'.^2);
diffd18Osd=diffd18Osd-diag(diag(diffd18Osd));

diffiv=diffd18O*scal;
diffivsd=sqrt((scalsd*diffd18O).^2+(diffd18Osd*scal).^2);

Csd=meshgrid(ivsd);
Csddiff=sqrt(Csd.^2+Csd'.^2); Csddiff=Csddiff-diag(diag(Csddiff));
diffivcv=.5*(Csddiff.^2-diffivsd.^2);
ivcv=diffivcv+diag(ivsd.^2);

iv0 = ivbase(end:-1:1); ivcv=ivcv(end:-1:1,end:-1:1);
t=LR04(end:-1:1,1);


T=meshgrid(t);
dT=T-T';
filt=exp(-(dT/timescale).^2);

%%%%%%

esl0=esl0(:);
icevol0=sum(esl0);

r0=mvnrnd(iv0,ivcv.*filt,Nsamp);

T=repmat(t(:)',[Nsamp 1]);

sub1=find(t>=ctrtime);
sub2=find(t<=ctrtime);
sub1=sub1(end:-1:1);

for j=1:Nsamp
	ctrshifted= ctrtime + ctrshiftsd*randn;
	dt=diff(t(sub1));
	keepdrawing=1;
	while keepdrawing
		shifts=randn(size(dt));
		dt2=dt+timestepsd*shifts;
		keepdrawing=sum((sign(dt(1))*dt2)<=.1);
		if (ctrshifted+sum(dt2))<max(t2)
			keepdrawing=1;
		end
	end
	T(j,sub1) = ctrshifted + [0 cumsum(dt2(:)')];

	dt=diff(t(sub2));
	keepdrawing=1;
	while keepdrawing
		shifts=randn(size(dt));
		dt2=dt+timestepsd*shifts;
		keepdrawing=sum((sign(dt(1))*dt2)<=.1);
		if (ctrshifted+sum(dt2))>min(t2)
			keepdrawing=1;
		end
	end
	dt=dt+timestepsd*shifts;
	T(j,sub2) = ctrshifted + [0 cumsum(dt2(:)')];
end

clear iced;
for j=1:Nsamp
	r(j,:)=interp1(T(j,:),r0(j,:),t2);
	sub=find(isnan(r(j,:)));
	if length(sub)>0
		disp('NaN found!'); pause
	end
	
	disp(['Sampling ' num2str(j) '/' num2str(Nsamp) '...']);
	iced(:,:,j)=SampleIceVolumePrior(icevol0+r(j,:),mu,cv,bins,t2);
end

% calculate gsl-varying bit of steric

% from Meehl et al. projections:
MeehlTE=[2 2 1.5 1.5 1.1 .9 .8 .6];
MeehlST=[4 4 4.5 2.5 3 4.5 2 3];
LGMSteric=(MeehlTE./MeehlST)*-5;

stericslop_mean = mean(LGMSteric/125);
stericslop_sd = std(LGMSteric/125);
stericslop_rnd = randn([Nsamp 1])*stericslop_sd + stericslop_mean;

% constant bit of steric
T2 = meshgrid(t2);
dT2 = T2-T2';
filt=exp(-(dT2/sterictimescale).^2);
stericcons=mvnrnd(zeros(size(t2)),stericsd*filt,Nsamp);

steric = repmat(stericslop_rnd,[1 size(r,2)]) .* r + stericcons;

gsl=-r+steric;



%%%%%%%
%%%%%%%
%%%%%%%

function LR04=LookupLR04(datset)


LR04_Holocene = [
0.0 3.23 0.03
1.0 3.23 0.04
2.0 3.18 0.03
3.0 3.29 0.03
4.0 3.30 0.03
5.0 3.26 0.03
6.0 3.33 0.04
7.0 3.37 0.04
8.0 3.42 0.03
9.0 3.38 0.04
10.0 3.52 0.04
11.0 3.60 0.04
12.0 3.92 0.05
13.0 4.06 0.04
14.0 4.28 0.03
15.0 4.49 0.04
16.0 4.75 0.03
17.0 4.88 0.04
18.0 5.02 0.03
19.0 4.96 0.03
20.0 4.99 0.04
21.0 4.91 0.03
22.0 4.88 0.03
23.0 4.86 0.03
24.0 4.81 0.04
25.0 4.82 0.02
26.0 4.67 0.04
27.0 4.75 0.04
28.0 4.75 0.03
];

LR04_LIG=[
80.0  3.90  0.04
81.0  3.82  0.04
82.0  3.80  0.04
83.0  3.83  0.04
84.0  3.82  0.04
85.0  3.95  0.06
86.0  4.06  0.05
87.0  4.18  0.05
88.0  4.11  0.05
89.0  4.08  0.04
90.0  4.06  0.05
91.0  4.03  0.04
92.0  3.98  0.04
93.0  3.90  0.03
94.0  3.84  0.04
95.0  3.77  0.05
96.0  3.75  0.03
97.0  3.83  0.04
98.0  3.90  0.04
99.0  3.78  0.06
100  3.81  0.06
101  3.92  0.05
102  3.86  0.05
103  3.88  0.04
104  3.92  0.05
105  3.85  0.06
106  4  0.04
107  4.04  0.05
108  4.11  0.04
109  4.12  0.05
110  4.04  0.05
111  4.02  0.05
112  4.03  0.04
113  3.93  0.05
114  3.81  0.05
115  3.71  0.04
116  3.58  0.04
117  3.54  0.05
118  3.44  0.06
119  3.3  0.03
120  3.27  0.04
121  3.26  0.04
122  3.18  0.03
123  3.1  0.05
124  3.27  0.04
125  3.14  0.06
126  3.16  0.04
127  3.37  0.05
128  3.71  0.07
129  3.9  0.08
130  3.67  0.09
131  3.81  0.09
132  4.2  0.06
133  4.41  0.05
134  4.7  0.05
135  4.86  0.04
136  4.82  0.05
137  4.8  0.05
138  4.89  0.04
139  4.87  0.05
140  4.98  0.05
141  4.81 0.05
142  4.75 0.06
143  4.78 0.05
144  4.82 0.05
145  4.74 0.05
146  4.77 0.05
147  4.82 0.05
148  4.71 0.05
149  4.75 0.04
150  4.75 0.06
151  4.66 0.06
152  4.64 0.05
153  4.62 0.04
154  4.66 0.04
155  4.51 0.04
156  4.78 0.04
157  4.74 0.04
158  4.69 0.03
159  4.66 0.05
160  4.69 0.05
];

LR04=eval(datset);