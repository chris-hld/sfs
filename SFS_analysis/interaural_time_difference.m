function [itd,idxleft,idxright] = interaural_time_difference(insigleft,insigright,fs,mode,fit)
%INTERAURAL_TIME_DIFFERENCE Extract the ITD between the two given signals
%   Usage: itd = interaural_time_difference(insigleft,insigright,fs)
%   	   itd = interaural_time_difference(insigleft,insigright,fs,mode)
%
%   Input parameters:
%       insigleft   - left ear signal. This can also be a matrix containing
%                     left signals for different frequency bands
%       insigright  - the same as insigleft, but for the right ear
%       fs          - sampling rate
%		mode		- describes which algorith was used: simple is default
%					- you should use 'int' too and compare results
%					- 'simple'		- fast, use only for low frequencies (default)
%					- 'hilbert' 	- use only for hight frequencies
%					- 'int'         - use if simple and hilbert fail
%					- 'int2'		- use for signals with high noise
%					- 'int3'		- use if int fails
%					- 'max'			- very primitive algorithm, doesn't work well
%		fit			- seconds the itd from angle i can difer from itd of angle i+1 
%		              (default is 0.1 )
%
%   Output parameters:
%       itd         - ITD for the given signals. A single value for two
%                     given signals or a vector with values for every
%                     frequency band
%
%   INTERAURAL_TIME_DIFFERENCE(insigleft,insigright,fs) extractes the ITD between
%   the left and right signal(s) by using an edge detection algorithm to identify
%   the first non-zero entry in both IRs and then calculating the time difference.
%
%R gaik1993, sandvad1994, lindau2010

% AUTHOR: Hagen Wierstorf,Lars-Erik Riechert


%% ===== Checking of input parameters ===================================
nargmin = 3;
nargmax = 5;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(insigleft,insigright);
isargpositivescalar(fs);
if size(insigleft)~=size(insigright)
    error('%s: insigleft and insigright have to be the same size!', ...
        upper(mfilename));
end
if nargin==5
    isargchar(mode);
elseif nargin==4
    isargchar(mode);
	fit=0.1;
elseif nargin==3
    mode='simple';
	fit=0.1;
end


%% ===== Computation ====================================================
if strcmp ( mode,'hilbert') 
	% Extract the envelope of the input signals
	insigleft = abs(hilbert(insigleft));
	insigright = abs(hilbert(insigright));
end

% See if we have more than one frequency channel in the insig
itd = zeros(1,size(insigleft,2));
%idxleft = itd;
%idxright = itd;
itdold =0;
jminold=0;
for ii = 1:size(insigleft,2)
if strcmp ( mode,'hilbert') |strcmp ( mode,'simple')
    % Treshold after sandvad1994 (5% of maximum in each IR)
    % NOTE: I have changed it to 10%
    tresholdleft = 0.15 * max(abs(insigleft(:,ii)));
    tresholdright = 0.15 * max(abs(insigright(:,ii)));
    % Ten fold upsampling (after lindau2010) to have a smoother output

    resampleft = resample(insigleft(:,ii),10*fs,fs);
    resampright = resample(insigright(:,ii),10*fs,fs);

	
	idxleft(ii) = find(abs(resampleft) > tresholdleft,1,'first');
    idxright(ii) = find(abs(resampright) > tresholdright,1,'first');
	itd(ii) = (idxleft(ii)-idxright(ii))/(10*fs);
	
	if 0 %mod(ii,30)==0
		range = 001:8000;
		plot(range,resampright(range),range,resampleft(range),idxleft(ii),tresholdleft,'+',idxright(ii),tresholdright,'*');
		grid on;
		file = sprintf('results/itd_%s_%i.png',mode,ii);
		print(file);
	end
	
elseif strcmp ( mode,'max')
    % Find maximum and use the maximum to calculate the ITD
    [maxleft idxleft(ii)] = max(abs(insigleft(:,ii)));
    [maxright idxright(ii)] = max(abs(insigright(:,ii)));
	itd(ii) = (idxleft(ii)-idxright(ii))/(10*fs);

elseif strncmp ( mode,'int',3)
	%here we try an other approach / integration of the difference and minimize
	difmin=0;
	
	a=sum(abs(insigleft(:,ii)))/sum(abs(insigright(:,ii)));			%here we normalize the both signals
	insigright(:,ii)=a*insigright(:,ii);							
																	%the right one should be louder
	%put zeros in the end and beginning
	l=(length(insigleft(:,ii)));
	%stretch= round(0.1*l);
	stretch = round(fs*0.001);
	
	newleft = zeros(1,l+2*stretch);
	newleft(stretch+1:l+stretch) = insigleft(:,ii);
	
	
	[maxright idxright] = max(abs(insigright(:,ii)));
	

	weight=(idxright*2-1:-1:idxright)-idxright;
	
	if strcmp(mode,'int2')	
		weight=log(weight+1);
	end
	
	%weight=(1-10*(abs((1:l+2*stretch)-jminold)/l));
	%for k = 1:l+2*stretch
	%if weight(k)<0
	%	weight(k)= 0;
	%end
	%end
	%weight=weight.^2;

	%add zeros also to the end and beginning of the right signal, and let the data turn from beginning to end
	
	for j = 1:2*stretch
		newright = zeros(1,l+2*stretch);						%create an empty signal of the new lenght
		if strcmp(mode,'int3')
			newright(j:j-1+l) = insigright(:,ii);
		else
			newright(j:j-1+idxright) = insigright(1:idxright,ii);	%we cut the right(shifting) signal after the maximum	
		end
		
		dif = newleft - newright;
		
		if strcmp(mode,'int') |strcmp(mode,'int2')
		    dif2 = zeros(1,l+2*stretch);
			dif2(j:j-1+idxright) = dif(j:j-1+idxright); 	%only the errors before we cut the signal are important
			%dif2(find(abs(newright)>0.05*maxright)) = dif(find(abs(newright)>0.05*maxright));	%only errors where the right signal is not zero are important
		
			dif=dif2;
			dif(j:j-1+idxright)=weight.*dif(j:j-1+idxright);	%here we say which errors are more important than others
		end
		
		%dif = dif.*weight*100;	
		%dif = dif.^2;
		dif=abs(dif);
		sumdif = sum(dif);


		if ((sumdif < difmin) & abs(((-stretch+j)/fs)-itdold)<fit   ) | j==1
			jmin = j;
			difmin=sumdif;
			difi=dif;
			itd(ii) = (-stretch+j)/fs;
		end
	end
	itdold=itd(ii);
	jminold=jmin;
	
	if 0 %mod(ii,30)==0																	%plots the signals for the given angels
		newright = zeros(1,l+2*stretch);
		if strcmp(mode,'int3')
			newright(jmin:jmin-1+l) = insigright(:,ii);
		else
			newright(jmin:jmin-1+idxright) = insigright(1:idxright,ii);	%we cut the right(shifting) signal after the maximum	
		end
		range = 200:1000;
		plot(range,newright(range),range,newleft(range),range,difi(range),jmin:jmin-1+idxright,weight*0.0000002);
		grid on;
		file = sprintf('results/itd_%s_%i.png',mode,ii);
		print(file);
	end

end
end
