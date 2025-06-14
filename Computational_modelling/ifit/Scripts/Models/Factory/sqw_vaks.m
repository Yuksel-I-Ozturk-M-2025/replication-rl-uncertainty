function signal=sqw_vaks(varargin)
% model = sqw_vaks(p, h,k,l,w, {signal}) : dispersion(HKL) in perovskites ABX3 with DHO(energy)
%
%   iFunc/sqw_vaks: a 4D S(q,w) with a HKL dispersion, and a DHO line
%      shape, based on a parameterisation by Vaks, suited for ABX3 perovskite 
%      compounds. This parameterisation is based on a quadratic expansion of the
%      Hamiltonian, valid for e.g. |k|<.3 rlu in a cubic 
%      crystal and provides energies of TA1,TA2,LA,TO1 and TO2 modes.
%
% WARNING:
%      Single intensity and line width parameters are used here.
%      The HKL position is relative to the closest Bragg peak. 
%
%   The acoustic parameters can be related to the elastic constants as:
%     At= (C12+2C44)/rho [meV2/rlu2]
%     Al=  C44/rho
%     Aa= (C11-C12-2C44)/rho
%   with rho=density of the crystal [g/cm3]
%
% Example:
% s=sqw_vaks('KTaO3'); qh=linspace(0,.5,50);qk=qh; ql=qh; w=linspace(0.01,10,51);
% f=iData(s,[],qh,qk,ql,w); scatter3(log(f(:,:,1,:)),'filled');
%
% References: E. Farhi et al, EPJB, 15 (2000) pp 615-623. DOI: 10.1007/s100510051164 
%             A.K. Tagantsev et al, Nature Comm, 4 (2013) 2229. DOI: 10.1038/ncomms3229
%             V.G. Vaks, Introduction to the Microscopic Theory of Ferroelectrics (Nauka, Moscow, 1973)
%
% input:  p: sqw_vaks model parameters (double)
%           p(1)=At transverse acoustic [meV2/rlu2]
%           p(2)=Al longitudinal acoustic [meV2/rlu2]
%           p(3)=Aa anisotropic acoustic [meV2/rlu2]
%           p(4)=St soft mode  transverse soft optical [meV2/rlu2]
%           p(5)=Sa soft mode  anisotropic soft optical [meV2/rlu2]
%           p(6)=Vt transverse acoustic-optical coupling [meV2/rlu2]
%           p(7)=Va anisotropic acoustic-optical coupling [meV2/rlu2]
%           p(8)=w0 soft mode frequency at q=0, depends on temperature [meV] 
%           p(9)=Gamma   dispersion DHO half-width in energy [meV]
%           p(10)=Temperature of the material [K]
%           p(11)=Amplitude
%           p(12)=Background (constant)
%          or 'KTaO3','BaTiO3','SrTiO3' for predefined settings
%         qh: axis along QH in rlu (row,double)
%         qk: axis along QK in rlu (column,double)
%         ql: axis along QL in rlu (page,double)
%         w:  axis along energy in meV (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
%
% Version: Aug. 22, 2017
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_sine3d
%   sqw_cubic_monoatomic, <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

signal.Name           = [ 'Sqw_Vaks phonon dispersion(HKL) for perovskites ABX3 [' mfilename ']' ];
signal.Description    = 'A phonon dispersion(HKL) for perovskites ABX3 with DHO(energy) shape. Parameterisation by Vaks.';

signal.Parameters     = {  ...
  'At transverse acoustic [meV2/rlu2]' ...
  'Al longitudinal acoustic [meV2/rlu2]' ...
  'Aa anisotropic acoustic [meV2/rlu2]' ...
  'St soft mode  transverse soft optical [meV2/rlu2]' ...
  'Sa soft mode  anisotropic soft optical [meV2/rlu2]' ...
  'Vt transverse acoustic-optical coupling [meV2/rlu2]' ...
  'Va anisotropic acoustic-optical coupling [meV2/rlu2]' ...
  'w0 soft mode frequency at q=0, depends on temperature [meV]' ...
  'Gamma Damped Harmonic Oscillator width in energy [meV]' ...
  'Temperature [K]' ...
  'Amplitude' 'Background' };
  
signal.Dimension      = 4;           % dimensionality of input space (axes) and result

signal.Guess          = [1553 4551.6 2264.9 4828.1 8956.1 2450.7 -3087.3 6.2 0.1 10 1 0];        % default parameters

% get code to read xyzt and build HKL list and convolve DHO line shapes
[script_hkl, script_dho] = sqw_phonons_templates;

% we create 5x5x(kx,ky,kz sets) matrices
signal.Expression     = { ...
  '% get parameter values',...
  'At=p(1); Al=p(2); Aa=p(3); St=p(4); Sa=p(5); Vt=p(6); Va=p(7); w0=p(8); L=w0^2;', ...
  script_hkl{:}, ...
  'HKL=HKL-round(HKL);', ....
  'FREQ=zeros(size(HKL,1),5);',...
  'for index=1:size(HKL,1); kx=HKL(index,1); ky=HKL(index,2); kz=HKL(index,3);', ...
  '  k2 = abs(kx^2 +  ky^2 +  kz^2); k = sqrt(k2);', ...
  '  n1=0; n2=0; n3=0; np=0; h11=0; h12=n1; h22=n1; h15=n1; h25=n1;', ...
  '  if k>1e-10, n1 = kx/k; n2 = ky/k; n3 = kz/k; else n1=0; n2=1; n3=0; end', ...
  '  np = sqrt(abs(1-n1*n1)); ',...
  '  if np>1e-10, h11 = 2*n2*n2*n3*n3/np/np;', ...
  '    h12 = -n1*n2*n3*(n3*n3-n2*n2)/np/np;', ...
  '    h22 = 2*n1*n1*(np*np-n2*n2*n3*n3/np/np);', ...
  '    h15 = -n2*n3*(n2*n2-n3*n3)/np;', ...
  '    h25 = n1*(n1*n1*np*np-n2^4-n3^4)/np;', ...
  '  else h11=0; h12=0; h22=0; h15=0; h25=0;',...
  '    n1=1; n2=0; n3=0; end',...
  '% now build the His and Hanis matrices', ...
  'His = [ (L+St*k2)		0		Vt*k2		0		0 ; 0		(L+St*k2)	0		Vt*k2		0 ; Vt*k2		0		At*k2		0		0 ; 0		Vt*k2		0		At*k2		0 ; 0		0		0		0		Al*k2 ];', ...
  'Hanis = eye(5);',...
  'Hanis(1,1) = Sa*h11;',...
	'Hanis(1,2) = Sa*h12;',...
	'Hanis(2,1) = Hanis(1,2);',...
	'Hanis(2,2) = Sa*h22;',...
	'Hanis(1,3) = Va*h11;',...
	'Hanis(1,4) = Va*h12;',...
	'Hanis(1,5) = Va*h15;',...
	'Hanis(2,3) = Va*h12;',...
	'Hanis(2,4) = Va*h22;',...
	'Hanis(2,5) = Va*h25;',...
	'Hanis(3,1) = Hanis(1,3);',...
	'Hanis(4,1) = Hanis(1,4);',...
	'Hanis(5,1) = Hanis(1,5);',...
	'Hanis(3,2) = Hanis(2,3);',...
	'Hanis(4,2) = Hanis(2,4);',...
	'Hanis(5,2) = Hanis(2,5);',...
	'Hanis(3,3) = Aa*h11;',...
	'Hanis(3,4) = Aa*h12;',...
	'Hanis(3,5) = Aa*h15;',...
	'Hanis(4,3) = Aa*h12;',...
	'Hanis(4,4) = Aa*h22;',...
	'Hanis(4,5) = Aa*h25;',...
	'Hanis(5,3) = Hanis(3,5);',...
	'Hanis(5,4) = Hanis(4,5);',...
	'Hanis(5,5) = Aa*(n1^4+n2^4+n3^4);	',...
  'Hanis = Hanis * k2; [eigvectors,eigvalues ] = eig(His + Hanis);',...
  '  eigvalues = sqrt(diag(eigvalues));',...
  '  [dummy,sorted]=max(abs(eigvectors''));',...
  '  eigvalues = eigvalues(sorted);',...
  '  FREQ(index,:) = eigvalues;', ...
  'end % index in kx,ky,kz', ...
  ' % multiply all frequencies(columns, meV) by a DHO/meV', ...
  'Gamma=p(9); T=p(10); Amplitude=p(11); Bkg=p(12); ', ...
  script_dho{:}, ...
 };

signal.UserData.properties.spacegroup        = 'Cubic (230)';
signal.UserData.properties.spacegroup_number = 230;

signal=iFunc(signal);

if nargin == 0
  varargin{1} = 'gui';
end
if numel(varargin) && ischar(varargin{1})
  switch lower(varargin{1})
  case 'gui'
    doc(iData,'Models.html#mozTocId226848');  % doc for Vaks
    NL = sprintf('\n');
    prompt = { [ '{\bf Enter Sqw Vaks crystal configuration}' NL ...
      'you can enter {\color{blue}KTaO3}, {\color{blue}SrTiO3} or {\color{blue}BaTiO3},' NL ...
      'or a {\color{blue}vector} of Vaks parameters (12 values) such as [1553 4551.6 2264.9 4828.1 8956.1 2450.7 -3087.3 6.2 0.1 10 1 0],' NL ...
      'or {\color{blue}defaults} to generate a KTaO3 Vaks model, ' NL ...
      'or {\color{blue}any expression} to evaluate and provide a 12 values-vector.' ]
    };
    dlg_title = 'iFit: Model: Sqw Vaks';
    defAns    = {'KTaO3'};
    num_lines = [ 1 ];
    op.Resize      = 'on';
    op.WindowStyle = 'normal';   
    op.Interpreter = 'tex';
    answer = inputdlg(prompt, dlg_title, num_lines, defAns, op);
    if isempty(answer), 
      signal = [];
      return; 
    end
    % now interpret the result
    answer = answer{1};
    NumEval = str2mat(answer);
    if isempty(NumEval)
      NumEval = answer;
      try
        if ~any(strcmpi(answer, {'KTaO3','SrTiO3','BaTiO3'}))
          NumEval = evalc(answer);
        end
      end
    end
    signal = sqw_vaks(NumEval);
    return
  case 'ktao3'
    St = 4828.100000; At = 1553.000000; Vt = 2450.700000;
	  Sa = 8956.100000; Aa = 2264.900000; Va = -3087.300000;
	  Al = 4551.600000;
	  disp([ mfilename ': using KTaO3 parameters' ])
	  p = [ At Al Aa St Sa Vt Va sqrt(6.2) .1 10 1 0 ];
	  signal.ParameterValues = p;
	  signal.Name = [ signal.Name ' KTaO3' ];
  case 'srtio3'
    St = 4*1080;  At = 1920; Vt = 2*1080; 
	  Sa = 25*1080; Aa = 1230; Va = -4860;
	  Al = 6500;
	  disp([ mfilename ': using SrTiO3 parameters' ])
	  p = [ At Al Aa St Sa Vt Va sqrt(4.2) .1 20 1 0 ];
	  signal.ParameterValues = p;
	  signal.Name = [ signal.Name ' SrTiO3' ];
  case 'batio3'
    St = 2*1080;  At = 1.8*1080; Vt = 0.1*1080; 
    Sa = 22*1080; Va = 5*1080;   Aa = -2.3*1080;
    Al = 5.8*1080; 
    disp([ mfilename ': using BaTiO3 parameters' ])
    p = [ At Al Aa St Sa Vt Va 1 .1 600 1 0 ];
    signal.ParameterValues = p;
    signal.Name = [ signal.Name ' BaTiO3' ];
  end
elseif nargin > 1
  signal = signal(varargin{:});
end

