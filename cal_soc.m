function soc = cal_soc(b, soc_init, yita_ch, yita_dc, B_E, ts)
N = length(b);
for i = 1:N+1
if(i==1)
    soc(1) = soc_init;
else
    if(b(i-1)>0) %last step discharge
        soc(i)=(soc(i-1)*B_E-(b(i-1)/yita_dc)*ts)/B_E;
    else % last step charge
        soc(i) = (soc(i-1)*B_E-b(i-1)*yita_ch*ts)/B_E;
    end
end
end