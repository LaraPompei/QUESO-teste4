cd outputData/
sipOutput_raw_chain
sipOutput_filt_chain

media_a = mean(ip_mh_filtChain_unified);
std_a = std(ip_mh_filtChain_unified);

disp("Media a recuperado: "), disp(media_a);
disp("Desvio Padrão a Recuperado: "), disp(std_a);

cd ../
a_data

media_a = mean(a_data);
std_a = std(a_data);

disp("Media a amostra: "), disp(media_a);
disp("Desvio Padrão a amostra: "), disp(std_a);
