function ppd = Secc_watson2018(S, ecc, pars)
%SECC_WATSON2018 Predicts ppd based on the Watson 2018 model.
%
%   ppd = SECC_WATSON2018(S, ecc, pars) computes the predicted pixels per
%   degree (ppd) using the sensitivity (S), eccentricity (ecc), and a set
%   of parameters (pars) based on the model proposed by Watson (2018):
%
%   Watson, A. B. (2018). The field of view, the field of resolution, and 
%   the field of contrast sensitivity. Electronic Imaging, 30, 1-11.
%
%   The function assumes a logarithmic relationship between sensitivity and
%   eccentricity, modulated by the parameters. The parameters pars(1),
%   pars(2), and pars(3) represent the offset log sensitivity, and scaling
%   factors cF and cR, respectively.
%
%   Inputs:
%       S - A vector of contrast sensitivity measurements.
%       ecc - A vector of eccentricity values corresponding to each
%             sensitivity measurement. Must be the same length as S.
%       pars - A vector of parameters [c, k1, k2] used in the model.
%
%   Outputs:
%       ppd - A vector of predicted pixels per degree, computed for each
%             pair of sensitivity and eccentricity values.


% assert(length(S) == length(ecc));

c = pars(1); % offset log sensitivity 
cF = pars(2);
cR = pars(3);

rho = (log10(S)-c)./(cF.*(1+(cR.*ecc)));
ppd = rho.*2;

end