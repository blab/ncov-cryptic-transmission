% random number distributions

function val = distribution(config, setSeedStr)

if nargin>1
    if ischar(config.simulation.seed) && strcmp(config.simulation.seed,'clock')
        rng(rem(now*100,1));
    else
        rng(config.simulation.seed);
    end
    return
end

switch config.distribution
    case 'fixed'
        val = config.value;
    case 'normal'
        val = config.params.mu + config.params.sigma*randn;
    case 'lognormal'
        val = exp(config.params.mu + config.params.sigma*randn);
    case 'positiveNormal'
        val = config.params.mu + config.params.sigma*randn;
        if val<0 
            val=0;
        end
    case 'negativeBinomial'
        p = 1/(1+config.params.mu/config.params.k);
        val = nbinrnd(config.params.k,p,1);
    case 'normalInteger'
        val = round(config.params.mu + config.params.sigma*randn);
    case 'lognormalInteger'
        val = ceil(exp(config.params.mu + config.params.sigma*randn));
    otherwise
        val = error('unknown distribution');
end

