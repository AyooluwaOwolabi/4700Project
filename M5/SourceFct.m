function E = SourceFct(t, InputParas)
    if isfield(InputParas, 'rep')
        n = floor(t / InputParas.rep);
        t = t - n * InputParas.rep;
    end

    if ~isstruct(InputParas)
        E = InputParas;
    else
        if isfield(InputParas, 'type') && strcmp(InputParas.type, 'square')
            % Square pulse: E0 for t0 < t < (t0 + wg), otherwise 0
            E = InputParas.E0 * (t >= InputParas.t0 & t <= (InputParas.t0 + InputParas.wg));
        else
            % Default: Gaussian pulse
            E = InputParas.E0 * exp(-(t - InputParas.t0)^2 / InputParas.wg^2) * ...
                exp(1i * (InputParas.we * t + InputParas.phi));
        end
    end
end
