function [ ipc,lifetime,energy,app_name ] = loaddata(  )

app_name = {'data_hmmer.txt' , 'data_mcf.txt', 'data_GemsFDTD.txt',...
         'data_lbm.txt'   , 'data_milc.txt', 'data_bwaves.txt',...
         'data_leslie3d.txt',  'data_stream.txt', 'data_gups.txt',...
         'data_libquantum.txt', 'data_zeusmp.txt'};
   
m  = length(app_name);


    for i = 1:m
        file_name = ['../',app_name{i}];
        fid = fopen(file_name);
        c = textscan(fid,'%s %f %f %f','HeaderLines', 1 ,'delimiter', ',');    
        fclose(fid);
    
        c{1,2} = c{1,2}/mean(c{1,2});
        c{1,3} = c{1,3}/mean(c{1,3});
        c{1,4} = c{1,4}/mean(c{1,4});
        ipc(:,i) = c{1,2}; 
        lifetime(:,i) = c{1,3}; 
        energy(:,i) = c{1,4}; 
        config = c{1,1};
    end

end

