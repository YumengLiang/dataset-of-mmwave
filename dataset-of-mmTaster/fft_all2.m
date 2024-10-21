%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: 
% Data extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author: Raghavendran
%%% Date: 24.04.2017
%%% Barbara Lenz on 04.07.2017: updated to new Matlab API
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [save_mat] = fft_all ( Rx1,Rx2,Rx3,Rx4,k1,k2,range)
%% Range calculation
Rx3 = transpose(Rx3);
Rx1 = transpose(Rx1);
Rx2 = transpose(Rx2);
Rx4 = transpose(Rx4);
mxFFT = [transpose(Rx3), transpose(Rx4) , transpose(Rx2) , transpose(Rx1)]; %64*4


mxMag = abs(mxFFT(k1:k2,:));
mxPhase = angle(mxFFT(k1:k2,:));

save_mat = zeros(2,length(k1:k2),4);
save_mat(1,:,:) = mxMag;
save_mat(2,:,:) = mxPhase;

end

