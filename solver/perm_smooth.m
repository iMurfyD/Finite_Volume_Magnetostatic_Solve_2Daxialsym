function p = perm_smooth(r,syst)
        s = smoothfn(r,syst);
        prat = syst.perm/syst.pfs;
        oneoverp = s*(1/prat) + (1-s)*1;
        p = 1/oneoverp;
        p = p*syst.pfs;
%         p = s*syst.perm + (1-s)*syst.pfs;
end
