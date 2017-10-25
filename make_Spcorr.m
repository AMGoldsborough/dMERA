function make_Spcorr(L,Jstr,Jdis,Jz,chi,Pdist,Jseed,ULmax,sweepmax,slow)
%make_Spcorr(L,Jstr,Jdis,Jz,chi,Pdist,Jseedmin,Jseedmax,ULmax,sweepmax)
% function to make spin corr from components
% input: Standard inputs
% output: files containing Spin correlation functions

% Andrew Goldsborough - 08/12/2016

% S.S = SzSz + 0.5*(SpSm + SmSp)

%open files to read in data
fnameSz = strcat('./Szcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_Szcorr_dMERA.txt');
Spcorr = importdata(fnameSz);

%check length
if size(Spcorr,1) == 0.5*L*(L-1)

    fnameSpSm = strcat('./SpSmcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_SpSmcorr_dMERA.txt');
    corr = importdata(fnameSpSm);
    Spcorr(:,3) = Spcorr(:,3) + 0.5*corr(:,3);
    
    fnameSmSp = strcat('./SmSpcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_SmSpcorr_dMERA.txt');
    corr = importdata(fnameSmSp);
    Spcorr(:,3) = Spcorr(:,3) + 0.5*corr(:,3);
    
    %open files to write to
    fnameSp = strcat('./Spcorr/',num2str(L),'_',num2str(Jstr),'_',num2str(Jdis),'_',num2str(Jz),'_',num2str(chi),'_',num2str(Pdist),'_',num2str(Jseed),'_',num2str(ULmax),'_',num2str(sweepmax),'_',num2str(slow),'_Spcorr_dMERA.txt');
    fprintf(strcat(fnameSp,'\n'));
    fidSpcorr = fopen(fnameSp, 'w');
    
    %print to file
    for i=1:size(Spcorr,1)
        fprintf(fidSpcorr,'%d %d %.15e\n',Spcorr(i,1),Spcorr(i,2),Spcorr(i,3));
    end
    
    %close file
    fclose(fidSpcorr);
else
    fprintf(fnameSz);
    fprintf(' : file not full\n');
end

