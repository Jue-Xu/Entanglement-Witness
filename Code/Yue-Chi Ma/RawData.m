classdef RawData < handle
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        raw_data_file
        raw_data_n = 0
        para
        vio_search
    end
    
    
    methods
        function rd = RawData(para)
            rd.para = para;
            rd.initialize();
        end
        
        function initialize(rd)
            rd.raw_data_file = rd.para.raw_data_file;
            para = rd.para;
            save(rd.raw_data_file,'para'); 
        end
        
        function data = gen_by_UPB(rd)
            p = rd.para;
            qubit_number = p.qubit_number;
            data_length = p.data_length;
            assert(qubit_number == 3,'qubit_number must be 3!');
            data = zeros(8,8,data_length);
            vecUPB = UPB([2,2,2]);
            rho = eye(8);
            for j = 1:4
                rho = rho - vecUPB(:,j)*vecUPB(:,j)';
            end    
            parfor i = 1:data_length
                A = Tensor(RandomInvertibleMatrix(2),...
                        RandomInvertibleMatrix(2),...
                        RandomInvertibleMatrix(2));
                rho_new = A*rho*A';
                data(:,:,i) = rho_new/trace(rho_new);
            end
            label = zeros(data_length,1);
            eval(['raw_data_',num2str(rd.raw_data_n),'= data;']);
            eval(['label_',num2str(rd.raw_data_n),'= label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)],'-append');
            rd.raw_data_n = rd.raw_data_n + 1;
        end
        
        function [data,label] = gen_by_Phi(rd)
            p = rd.para;
            disp('gen_multi_etg...');
            [data,label] = rd.gen_data_by_phi_rotation;
            disp('save raw data with different Phi...');
            raw_data_0 = data;
            label_0 = label;
            eval(['raw_data_',num2str(rd.raw_data_n),'= data;']);
            eval(['label_',num2str(rd.raw_data_n),'= label;']);
            save(rd.raw_data_file,'raw_data_0','label_0','-append');
        end
        
        function data = gen_by_BiSep(rd)
            p = rd.para;
            qubit_number = p.qubit_number;
            data_length = p.data_length;
            assert(qubit_number == 3,'qubit_number must be 3!');
            for mode = 0:2
                [data,label] = rd.gen_qu3_bisep(data_length,mode);
              %  label = ones(data_length,1)*mode;
                label = label*(mode + 1);
                disp(['sum of ent/all = (0.75 best) ',num2str(sum(label)/length(label)/(mode + 1))])
                disp(['save qu3 sep, mode = ',num2str(mode)]);
                eval(['raw_data_',num2str(mode),'= data;']);
                eval(['label_',num2str(mode),'= label;']);
                save(rd.raw_data_file,['raw_data_',num2str(mode)],...
                    ['label_',num2str(mode)],'-append');
            end
            
            
        %    rd.raw_data_n = rd.raw_data_n + 3;
        end
        
        function [data,neg] = gen_by_PureState(rd)
            if ~exist('label_for3','var')
                label_for3 = rd.raw_data_n;
            end
            p = rd.para;
            disp('gen_qu_sep...');
            data = rd.gen_qu_depolarize_etg(0);   
            data_length = p.data_length;
            label = ones(data_length,1)*label_for3;
            disp('save qu3 sep...');
            eval(['raw_data_',num2str(rd.raw_data_n),'= data;']);
            eval(['label_',num2str(rd.raw_data_n),'= label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)],...
                ['neg_',num2str(rd.raw_data_n)], '-append');
       
            rd.raw_data_n = rd.raw_data_n + 1;
        end
 
        function gen_by_Uni_Rotation_and_PPT_for_3(rd)
            p = rd.para;
            qubit_number = p.qubit_number;
            data_length = p.data_length;
            assert(qubit_number == 3,'qubit_number must be 3!');
            [data,label] = rd.gen_by_unitary_rotation_detect(data_length);
            eval(['raw_data_',num2str(rd.raw_data_n),'= data;']);
            eval(['label_',num2str(rd.raw_data_n),'= label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)],'-append');
            rd.raw_data_n = rd.raw_data_n + 1;
        end
        
        function data = gen_by_EntanglementWithNoise(rd)
            p = rd.para;
  %          qubit_number = p.qubit_number;
            data_length = p.data_length;
%            assert(qubit_number == 3,'qubit_number must be 3!');
            data = rd.gen_qu_etg_with_sep_noise();
            label = ones(data_length,1)*rd.raw_data_n;
            eval(['raw_data_',num2str(rd.raw_data_n),'= data;']);
            eval(['label_',num2str(rd.raw_data_n),'= label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)],'-append');
            rd.raw_data_n = rd.raw_data_n + 1;
        end 
        
        function data = gen_by_Division(rd,label_for3)
            if ~exist('label_for3','var')
                label_for3 = rd.raw_data_n;
            end
            p = rd.para;
 %           qubit_number = p.qubit_number;
 %           assert(qubit_number == 3,'qubit_number =3!');
            disp('gen_qu_sep...');
            data = rd.gen_qu_sep();
            data_length = p.data_length;
            label = ones(data_length,1)*label_for3;
            disp('save qu3 sep...');
            eval(['raw_data_',num2str(rd.raw_data_n),'= data;']);
            eval(['label_',num2str(rd.raw_data_n),'= label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)],'-append');
            rd.raw_data_n = rd.raw_data_n + 1;
        end 
       
        function data = gen_pure_by_locc(rd,data_length, pure_state)
            D = length(pure_state);
            qubit_number = rd.para.qubit_number;
      %      data_fully_sep = rd.gen_qu_sep(data_length);
            data = nan(D,D,data_length);

            parfor i = 1:data_length
                rho = pure_state*pure_state';
                A = RandomUnitary(2); %RandomInvertibleMatrix(2);
                for j = 1:qubit_number-1
                    A = kron(A, RandomUnitary(2));%RandomInvertibleMatrix(2));
                end
                rho_new = A*rho*A';
                data_ent = rho_new/trace(rho_new);
              %  p = rand();
                data(:,:,i) = data_ent; % + (1 - p) * data_fully_sep(:,:,i);
            end
        end
        
        function gen_GHZ(rd)
            p = rd.para;
            data_length = p.data_length;
            z = zeros(2^p.qubit_number,1);
            z(1) = 1;
            z(end) = 1;
            pure_state = z/sqrt(2);
            data = rd.gen_pure_by_locc(data_length, pure_state);
            label = zeros(data_length,1);
            eval(['raw_data_',num2str(rd.raw_data_n),'= data;']);
            eval(['label_',num2str(rd.raw_data_n),'= label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)]);
            rd.raw_data_n = rd.raw_data_n + 1;
        end
        
        function [data,label] = gen_data_by_phi_rotation(rd)
            p = rd.para;
            qubit_number = p.qubit_number;
            data_length = p.data_length;
            D = 2^qubit_number;
            data = zeros(D,D,data_length);
            label = -ones(data_length,1);
%             switch qubit_number
%                 case 2
%                 otherwise
%                     error('qubit_number must be 2 !!');
%             end
            max_mix = eye(D)/D;
            % sys = ones(1,qubit_number)*2;
            
           assert(qubit_number<3.5,'if qubit_number == 4, I should consider the p_c strategy!!')
           parfor i = 1:data_length
                th = rand()*pi/2;
                phi = rand()*2*pi; 
                state = [0,cos(th),exp(phi*1i)*sin(th),0]';
                rho = state * state';
                p = mod(randn()*0.02 + 1/(1 + D * cos(th)*sin(th)),1);  % randn()*0.02 ！！ 
                data(:,:,i) = p*rho + (1-p) * max_mix;  % Be careful!!
                label(i) = IsPPT(data(:,:,i));
           end
        end
                 
        function etg_data = gen_qu_depolarize_etg(rd)
            disp('gen_qu_depolarize_etg')
            p = rd.para;
            D = 2^p.qubit_number;
            assert D == 16
            data_length = p.data_length;
            etg_data = zeros(D,D,data_length);
     %       neg = nan(data_length,1);
            parfor nn=1:data_length
                state = RandomStateVector(D);
                etg_data(:,:,nn) = state * state';
%                 n1 = min(eig(PartialTranspose(etg_data(:,:,nn),1,[2,2,2,2])));
%                 n2 = min(eig(PartialTranspose(etg_data(:,:,nn),2,[2,2,2,2])));
%                 n3 = min(eig(PartialTranspose(etg_data(:,:,nn),3,[2,2,2,2])));
%                 n4 = min(eig(PartialTranspose(etg_data(:,:,nn),4,[2,2,2,2])));
%                 neg(nn)=min(n1,n2,n3,n4);
            end
        end
        
        function etg_data = gen_qu_etg_with_sep_noise(rd)
            p = rd.para;
            max_noise = p.max_noise;
            D = 2^p.qubit_number;
            data_length = p.data_length;
            disp('first generate sep data');
            tic
            sep_data = rd.gen_qu_sep();
            toc
            disp(['then generate ent data by add sep max noise = ',num2str(max_noise)]);
            etg_data = zeros(D,D,data_length);
            parfor nn = 1:data_length
                noi = rand()*max_noise;
                state = RandomStateVector(D);
                etg_data(:,:,nn) = (1 - noi) * (state * state') + noi * sep_data(:,:,nn);
            end
        end
        
        function sep_data = gen_qu_sep_OLD(rd)  % a little slow
            p = rd.para;
            d = p.qubit_number;
            D=2^d;
            data_length = p.data_length;
            sep_data = zeros(D,D,data_length);
            div = p.div;      
            %parfor
            for n = 1:num
                p = rand(div,1);
                p = p./sum(p);
                for dis = 1:div
                        den = kron(RandomDensityMatrix(2),...
                            RandomDensityMatrix(2));
                        for di = 1:(d-2)
                            den = kron(den,RandomDensityMatrix(2));
                        end
                    sep_data(:,:,n) = sep_data(:,:,n) + p(dis) * den;
                end
            end
        end
        
        function sep_data = gen_qu_sep(rd,data_length)
            p = rd.para;
            qnum = p.qubit_number;
            D = 2^qnum;
            if ~exist('data_length','var')
                data_length = p.data_length;
            end
            div = p.div;
            
            one_pure_data = ...
                rd.gen_data_by_sep_random_matrix(data_length*div,qnum); 
            sep_data = reshape(one_pure_data,[D,D,data_length,div]);
            p = rand(data_length,div);
            p = p./repmat(sum(p,2),1,div);
            p = reshape(p,[1,1,data_length,div]);
            disp('now bsxfun will work');
            sep_data = bsxfun(@times,p,sep_data); %mimic kron
            sep_data = sum(sep_data,4);    
        end  %faster than gen_qu_sep_OLD
        
        function data = gen_full_sep_for_UPB_OLD(rd) % for UPB, may be old in future
            p = rd.para;
            disp('gen_qu_sep...');
            data = rd.gen_qu_sep();
            data_length = p.data_length;
            label = ones(data_length,1);
            eval(['raw_data_',num2str(rd.raw_data_n),'= data;']);
            eval(['label_',num2str(rd.raw_data_n),'= label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)],'-append');
            rd.raw_data_n = rd.raw_data_n + 1;
        end
        
        function gen_Bisep_for_UPB(rd) % For UPB
            disp('generate biseparable data with multi-mode')
            pa = rd.para;
            data_length = pa.data_length;
            bisep_data_1 = rd.gen_general_qu3_bisep_for_UPB(data_length,1);
            bisep_data_2 = rd.gen_general_qu3_bisep_for_UPB(data_length,2);
            bisep_data_3 = rd.gen_general_qu3_bisep_for_UPB(data_length,3);
            bisep_data_mul = nan(8,8,data_length);
            disp('for par loop ...')
            parfor i = 1:data_length
                p = rand(3,1);
                p = p./sum(p);
                bisep_data_mul(:,:,i) = p(1)*bisep_data_1(:,:,i) + ...
                    p(2)*bisep_data_2(:,:,i) + ...
                    p(3)*bisep_data_3(:,:,i);
            end
            label = ones(data_length,1);
            eval(['raw_data_',num2str(rd.raw_data_n),'= bisep_data_mul;']);
            eval(['label_',num2str(rd.raw_data_n),'= label;']);
            save(rd.raw_data_file,['raw_data_',num2str(rd.raw_data_n)],...
                ['label_',num2str(rd.raw_data_n)],'-append');
            rd.raw_data_n = rd.raw_data_n + 1;
        end
        
        function bisep_data = gen_general_qu3_bisep_for_UPB(~,num,mode)
            disp(['generate general biseparable data (for UPB) with mode = ',num2str(mode)])
            d = 3;
            D = 2^d;
            bisep_data = zeros(D,D,num);
            parfor n = 1:num
                 bisep_data(:,:,n) = kron(RandomDensityMatrix(4),RandomDensityMatrix(2));
                 switch mode
                     case 0
                         bisep_data(:,:,n)=PermuteSystems(bisep_data(:,:,n),[3,1,2]);
                     case 1
                         bisep_data(:,:,n)=PermuteSystems(bisep_data(:,:,n),[1,3,2]);
                     case 2
                 end
            end
        end           
        
        function [bisep_data,label] = gen_qu3_bisep(~,num,mode)
            disp(['generate biseparable data with mode = ',num2str(mode)])
            d = 3;
            D = 2^d;
            bisep_data = zeros(D,D,num);
            div = 17;
            label = zeros(num,1);
        %    max_noise = 0.76;
            parfor n = 1:num
                 p = rand(div,1);
                 p = p./sum(p);
                 sep_den = zeros(8,8);
                 for dis = 1:div
                     sep_den = sep_den + p(dis)*kron(RandomDensityMatrix(2),kron(RandomDensityMatrix(2),RandomDensityMatrix(2)));
                 end
                 state = RandomStateVector(4);
                 noi = rand();
                 bisep_data(:,:,n) = (1-noi)*kron(state * state',RandomDensityMatrix(2)) + noi*sep_den; %iden to RDM(2)
                 label(n) = 1-IsPPT(PartialTrace(bisep_data(:,:,n),3,[2,2,2]));
                 switch mode
                     case 0
                         bisep_data(:,:,n)=PermuteSystems(bisep_data(:,:,n),[3,1,2]);
                     case 1
                         bisep_data(:,:,n)=PermuteSystems(bisep_data(:,:,n),[1,3,2]);
                     case 2
                 end
            end
        end   
        
        function [bisep_data,label] = gen_qu3_bisep_OLD(~,num,mode)
            disp(['generate biseparable data with mode = ',num2str(mode)])
            d = 3;
            D = 2^d;
            bisep_data = zeros(D,D,num);
            label = zeros(num,1);

          %  ini_state = [1,0,0,1]'/sqrt(2);
            div = 3;
         %   noise = 0;
            max_mix = eye(4)/4;
            %parfor
            for n = 1:num
                p = rand(div,1);
                p = p./sum(p);
                    for dis = 1:div
                       % state = kron(RandomUnitary(2),RandomUnitary(2))*ini_state;
                        state = RandomStateVector(4);
                    %    noi = rand(1)*noise;
                    %    ent = (1 - noi) * (state * state') + noi * max_mix;
                        noise = rand()*0.3;
                        ent = (1 - noise) * (state * state') + noise * max_mix;
                        bisep_data(:,:,n) = bisep_data(:,:,n) + p(dis) * kron(ent,RandomDensityMatrix(2));
                    end
                    label(n) = 1-IsPPT(PartialTrace(bisep_data(:,:,n),3,[2,2,2]));
                 switch mode
                     case 0
                         bisep_data(:,:,n)=PermuteSystems(bisep_data(:,:,n),[3,1,2]);
                     case 1
                         bisep_data(:,:,n)=PermuteSystems(bisep_data(:,:,n),[1,3,2]);
                     case 2
                 end
            end
        end   %old method
                
        function data = gen_data_by_sep_random_matrix(~,data_length,d) %old method?
            assert(d > 1.5,'d must >= 2 !!!');
            D = 2^d;
            data = zeros(D,D,data_length);
            parfor nn = 1:data_length
                den = kron(RandomDensityMatrix(2),...
                    RandomDensityMatrix(2));
                for di = 1:(d-2)
                    den = kron(den,RandomDensityMatrix(2));
                end
                data(:,:,nn) = den;
            end
        end    
        
        
    end
    
    methods(Static)
        function PPT_detection = check_by_PPT(raw_data)
            data_size = size(raw_data);
            D = data_size(1);
            data_length = data_size(3);
            PPT_detection = zeros(data_length,1);
            d = log2(D);
            sys = ones(1,d)*2;
            switch d
                case 4
                    parfor i = 1:data_length
                        PPT_detection(i) = ...
                            min(real(eig(PartialTranspose(raw_data(:,:,i),1,sys))))<0||...
                            min(real(eig(PartialTranspose(raw_data(:,:,i),2,sys))))<0||...
                            min(real(eig(PartialTranspose(raw_data(:,:,i),3,sys))))<0||...
                            min(real(eig(PartialTranspose(raw_data(:,:,i),4,sys))))<0
                        PPT_detection(i) = 1 - PPT_detection(i);
                    end
                case 3
                    parfor i = 1:data_length
                        PPT_detection(i) = ...
                            min(real(eig(PartialTranspose(raw_data(:,:,i),1,sys))))<0||...
                            min(real(eig(PartialTranspose(raw_data(:,:,i),2,sys))))<0||...
                            min(real(eig(PartialTranspose(raw_data(:,:,i),3,sys))))<0
                        PPT_detection(i) = 1 - PPT_detection(i);
                    end
            end
        end
        
        function PPT_detection = check_by_PPT_BAD(raw_data)
            data_size = size(raw_data);
            D = data_size(1);
            data_length = data_size(3);
            PPT_detection = zeros(data_length,1);
            d = log2(D);
            sys = ones(1,d)*2;
            switch d
                case 4
                    parfor i = 1:data_length
                        PPT_detection(i) = ...
                            min(eig(PartialTranspose(raw_data(:,:,i),1,sys)))<0||...
                            min(eig(PartialTranspose(raw_data(:,:,i),2,sys)))<0||...
                            min(eig(PartialTranspose(raw_data(:,:,i),3,sys)))<0||...
                            min(eig(PartialTranspose(raw_data(:,:,i),4,sys)))<0
                        PPT_detection(i) = 1 - PPT_detection(i);
                    end
                case 3
                    parfor i = 1:data_length
                        PPT_detection(i) = ...
                            min(eig(PartialTranspose(raw_data(:,:,i),1,sys)))<0||...
                            min(eig(PartialTranspose(raw_data(:,:,i),2,sys)))<0||...
                            min(eig(PartialTranspose(raw_data(:,:,i),3,sys)))<0
                        PPT_detection(i) = 1 - PPT_detection(i);
                    end
            end
        end
    end
end

