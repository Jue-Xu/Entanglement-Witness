classdef GapRawData < RawData
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        negativity
        gap_raw_data
        label
    end
    
    methods
        function grd = GapRawData(para)
            grd = grd@RawData(para);
        end
        
        function initialize(rd)
            para = rd.para;
            rd.raw_data_file = para.raw_data_file;    
            save(rd.raw_data_file,'para');
        end
        
        function gen_by_Gap(rd)
            p=rd.para;
        %    neg = p.negativity;
            disp('gen_data_by_Gap...')
            rd.gen_data();
            sep_in_all = sum(rd.negativity>0)/length(rd.negativity);
            negativity = rd.negativity;
            disp(['the ratio of separable and all data is: ',num2str(sep_in_all)]);

            eval(['raw_data_',num2str(rd.raw_data_n),'= rd.gap_raw_data;']);
            eval(['label_',num2str(rd.raw_data_n),'= rd.label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)],...
                'negativity','-append');
            rd.raw_data_n = rd.raw_data_n + 1;            
        end
        
        function gen_by_Diri(rd)
            p=rd.para;
        %    neg = p.negativity;
            disp('gen_data_by_Diri...')
            rd.gen_diri_data();
            if p.fill == true
                rd.filter(p.min_neg,p.max_neg) % neg min = -0.5, max = 2.5
            end
            sep_in_all = sum(rd.negativity>0)/length(rd.negativity);
            negativity = rd.negativity;
            disp(['the ratio of separable and all data is: ',num2str(sep_in_all)]);

            eval(['raw_data_',num2str(rd.raw_data_n),'= rd.gap_raw_data;']);
            eval(['label_',num2str(rd.raw_data_n),'= rd.label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)],...
                'negativity','-append');
            rd.raw_data_n = rd.raw_data_n + 1;            
        end
        
        function gen_data(rd)
            p = rd.para;
            num = p.data_length;
            qnum = p.qubit_number;
            
            
            
            tic  %parfor or for
            switch qnum
                case 2
                    D = 2^qnum;
                    a = ones(1,D)*0.5;
                    my_data = nan(D,D,num);
                    negativity = nan(num,1);
                    parfor nn = 1:num
                        mt = RandomDensityMatrix(D);
                        % U = RandomUnitary(D);
                        % mt = U*diag(drchrnd(a,1))*U'
                        my_data(:,:,nn) = mt;
                       % [my_label(nn),witness] = IsPPT(mt);
                       % rd.wit(nn) = real(witness'*PartialTranspose(mt)*witness);
                        negativity(nn) = min(real(eig(PartialTranspose(mt))));
                    end
                case 2.5
                    D = 6;
                    a = ones(1,D)*0.5;
                    sys_dim = [2,3];
                    for nn = 1:num
                       U = RandomUnitary(D);
                       mt = U*diag(drchrnd(a,1))*U';
                       my_data(:,:,nn) = mt;
                       neg_min = inf;
                       for q = 1:qnum
                           neg_min = min(neg_min,min(real(eig(PartialTranspose(mt,q,sys_dim)))));    
                       end
                       negativity(nn) = neg_min;
                    end
                case 3
                    D = 8;
                    my_data = nan(D,D,num);
                    negativity = nan(num,1);
                    sys_dim = [2,2,2];
                    a = ones(1,D)*0.5;
                    for nn = 1:num
                       U = RandomUnitary(D);
                       mt = U*diag(drchrnd(a,1))*U';
                       my_data(:,:,nn) = mt;
                       neg_min = inf;
                       for q = 1:3
                           neg_min = min(neg_min,min(real(eig(PartialTranspose(mt,q,sys_dim)))));    
                       end
                       negativity(nn) = neg_min;
                     %  negativity(nn) = min(eig(PartialTranspose(mt,1,[2,2,2])));
                    end
            end
            toc
            
            rd.negativity = negativity ;
            rd.label = negativity>=0;
            % check = my_label ~= rd.label;
            % sum(check)
            rd.gap_raw_data = my_data;
        end
        
        function gen_diri_data(rd)
            p = rd.para;
            num = p.data_length;
            qnum = p.qubit_number;
            D = 2^qnum;
            my_data = nan(D,D,num);
            negativity = nan(num,1);
            a = ones(1,D)*0.5;
         %   r = drchrnd(a,num);
            tic  %parfor or for
            
            parfor nn = 1:num
                U = RandomUnitary(D);
                mt = U*diag(drchrnd(a,1))*U';
                my_data(:,:,nn) = mt;
                neg_min = inf;
                sys_dim = ones(1,qnum)*2;
                for q = 1:qnum
                    neg_min = min(neg_min,min(eig(PartialTranspose(mt,q,sys_dim))));    
                end
                negativity(nn) = neg_min;
            end
            toc
            
            rd.negativity = negativity;
            rd.label = negativity>=0;
            rd.gap_raw_data = my_data;
        end
    end
end

