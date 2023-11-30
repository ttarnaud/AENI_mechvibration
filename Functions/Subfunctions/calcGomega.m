function Out = calcGomega(freqs,iSum,RScalp,f1,f2,sigmaSkull,sigmaGM,epsilonSkull,epsilonGM,norm_PosCDInit)
eps0 = 8.85*10^-12;
if isrow(freqs)
    freqs = freqs';
elseif ~iscolumn(freqs)
    error('freqs')
end    
CsigmaSkull = complex(sigmaSkull(abs(freqs)),2*pi*freqs.*epsilonSkull(abs(freqs)).*eps0);
CsigmaGM = complex(sigmaGM(abs(freqs)),2*pi*freqs.*epsilonGM(abs(freqs)).*eps0);
RatioSkT = CsigmaSkull./CsigmaGM;
gi = ((iSum+1).*RatioSkT+iSum).*(iSum.*RatioSkT./(iSum+1)+1)+(1-RatioSkT).*((iSum+1).*RatioSkT+iSum).*(f1.^(2*iSum+1)-f2.^(2*iSum+1))...
    -iSum.*(1-RatioSkT).^2.*(f1/f2).^(2*iSum+1);
ABC = RatioSkT*(2*iSum+1).^3./(gi.*(iSum+1).*iSum);
Out = 1./(4*pi.*CsigmaGM.*RScalp.^2).*ABC.*(norm_PosCDInit./RScalp).^(iSum-1);
Out = Out';
end