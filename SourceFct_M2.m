function E = SourceFct_M2(t,InputParas)

if isfield(InputParas, 'rep')
    n = floor(t/InputParas.rep);
    t = t-n*InputParas.rep;
end 

if ~isstruct(InputParas)
    E = InputParas;
else 
    E = InputParas.E0*exp(-(t-InputParas.t0)^2/InputParas.wg^2)*exp(1i*(InputParas.we*t + InputParas.phi));

end 

end 
