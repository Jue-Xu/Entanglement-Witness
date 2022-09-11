classdef QData < handle
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        para
        fea_num
        mat_dir = 'mat/'
        raw_data
        label
        current_label
        raw_data_file
        data_file
        data
        witness_list
        negativity
        si = eye(2);
        sx=[0,1;1,0];
        sy=[0,-1i;1i,0];
        sz=[1,0;0,-1];
    end
    
    methods
        function q = QData(para)
            q.para = para;
        end
        
        function parse_raw_data(q,do_load_neg)
            q.raw_data_file = q.para.raw_data_file;
            load(q.raw_data_file);
            
            q.raw_data = raw_data_0;
            q.label = label_0;
            
            index = 1;
            while(exist(strcat('raw_data_' , num2str(index)),'var'))
                eval(['q.raw_data = cat(3,q.raw_data, raw_data_',num2str(index),');']);
                eval(['q.label = [q.label; label_',num2str(index),'];']);
                index = index + 1;
            end
            
            if exist('do_load_neg','var') && do_load_neg == 1
                q.negativity = negativity;
            end
        end
       
        function witness_list = gen_CHSH_wit(q,theta)
             a = q.sz;
             if ~exist('theta','var')
                 theta = pi/4;
             end
             b = -cos(theta) * q.sz + sin(theta) * q.sx;
             aa = q.sx;
             bb = cos(theta) * q.sz + sin(theta) * q.sx;
%             a = (q.sx+q.sz)/sqrt(2);
%             aa = (q.sy+q.sz)/sqrt(2);
%             b = (q.sx-q.sz)/sqrt(2);
%             bb = (q.sy-q.sz)/sqrt(2);
            witness_list = zeros(4,4,4);
            witness_list(:,:,1) = kron(a,b);
            witness_list(:,:,2) = kron(a,bb);
            witness_list(:,:,3) = kron(aa,b);
            witness_list(:,:,4) = kron(aa,bb);
            q.witness_list = witness_list;
            q.data_file = [q.mat_dir , ['CHSH_like_data_for_2_qubits']];
        end
        
        function witness_list = gen_ShangHaiJiaoTong_wit(q)
             a = [sqrt(3)/2,1/2;1/2,-sqrt(3)/2];
             b = [1/2,sqrt(3)/2;sqrt(3)/2,-1/2];
             aa = (q.sy + q.sz)/sqrt(2);
             bb = (q.sy - q.sz)/sqrt(2);
            witness_list = zeros(4,4,4);
            witness_list(:,:,1) = kron(a,b);
            witness_list(:,:,2) = kron(a,bb);
            witness_list(:,:,3) = kron(aa,b);
            witness_list(:,:,4) = kron(aa,bb);
            q.witness_list = witness_list;
            q.data_file = [q.mat_dir , 'CHSH_like_data_for_2_qubits'];
        end   
        
        function val = get_CHSH_inequality_value(~,pure_state)
            den = pure_state*pure_state';
            v = zeros(1,4);
            for j=1:4
                v(j) = trace(den*q.witness_list(:,:,j)); 
            end
            val = sum(v);
        end
        
        function witness_list = gen_Tomography_wit(q)
            singlePauli = zeros(2,2,4);
            singlePauli(:,:,1)=q.si;
            singlePauli(:,:,2)=q.sx;
            singlePauli(:,:,3)=q.sy;
            singlePauli(:,:,4)=q.sz;
            
            qnum = q.para.qubit_number;
            D = 2^qnum;
            if qnum == 2.5
                D = 6;
            end
            witness_list = zeros(D,D,4^qnum);
            if qnum == 2
                for ii =1:4
                    for jj=1:4
                        witness_list(:,:,jj+(ii-1)*4)=kron(singlePauli(:,:,ii),singlePauli(:,:,jj));
                    end
                end
            elseif qnum == 3
                for ii =1:4
                    wit1 = singlePauli(:,:,ii);
                    for jj=1:4
                        wit2 = kron(wit1,singlePauli(:,:,jj));
                        for kk=1:4
                            witness_list(:,:,(ii-1)*16+(jj-1)*4+kk) = kron(wit2,singlePauli(:,:,kk));
                        end
                    end
                end
            elseif qnum == 4
                for ii =1:4
                    wit1 = singlePauli(:,:,ii);
                    for jj=1:4
                        wit2 = kron(wit1,singlePauli(:,:,jj));
                        for kk=1:4
                            wit3 = kron(wit2,singlePauli(:,:,kk));
                            for ll = 1:4
                                witness_list(:,:,(ii-1)*64+(jj-1)*16+(kk-1)*4+ll) = kron(wit3,singlePauli(:,:,ll));
                            end
                        end
                    end
                end
            else
                error('qubit_number must be 2, 3 or 4!!!');
            end
                
            q.witness_list = witness_list(:,:,2:end);
            q.data_file = [q.mat_dir , ['FTI_data_for_',num2str(qnum),'_qubits']];
        end 
        
        function witness_list = gen_CHSH_like_wit(q,fea_group_num,fea_num)
            if ~exist('fea_num','var')
                fea_num = 4;
            end
            q.fea_num = fea_num;
            qnum = q.para.qubit_number;
            disp(['generate ',num2str(fea_group_num),' groups of CHSH-like features']);
            
            witness_list = zeros(2^qnum,2^qnum,fea_group_num * fea_num);
            switch qnum
                case 2
                    switch fea_num
                        case 4
                            for ii = 1:fea_group_num
                                a = q.gen_local_random_observe();
                                b = q.gen_local_random_observe();
                                aa = q.gen_local_random_observe();
                                bb = q.gen_local_random_observe();

                                witness_list(:,:,ii) = kron(a,b);
                                witness_list(:,:,ii + fea_group_num) = kron(aa,bb);
                                witness_list(:,:,ii + fea_group_num*2) = kron(a,bb);
                                witness_list(:,:,ii + fea_group_num*3) = kron(aa,b);
                            end  
                        case 8
                            for ii = 1:fea_group_num
                                a = q.gen_local_random_observe();
                                b = q.gen_local_random_observe();
                                aa = q.gen_local_random_observe();
                                bb = q.gen_local_random_observe();

                                witness_list(:,:,ii) = kron(a,b);
                                witness_list(:,:,ii + fea_group_num) = kron(aa,bb);
                                witness_list(:,:,ii + fea_group_num*2) = kron(a,bb);
                                witness_list(:,:,ii + fea_group_num*3) = kron(aa,b);

                                witness_list(:,:,ii + fea_group_num*4) = kron(aa,q.si);
                                witness_list(:,:,ii + fea_group_num*5) = kron(a,q.si);
                                witness_list(:,:,ii + fea_group_num*6) = kron(q.si,b);
                                witness_list(:,:,ii + fea_group_num*7) = kron(q.si,bb);
                            end  
                        otherwise
                            error('fea_num must be 4 or 8 !!!');
                    end
                case 3
                    switch fea_num
                        case 4
                            for ii = 1:fea_group_num
                                a = q.gen_local_random_observe();
                                b = q.gen_local_random_observe();
                                c = q.gen_local_random_observe();
                                aa = q.gen_local_random_observe();
                                bb = q.gen_local_random_observe();
                                cc = q.gen_local_random_observe();

                                witness_list(:,:,ii) = kron(kron(a,b),c);
                                witness_list(:,:,ii + fea_group_num) = kron(kron(aa,bb),c);
                                witness_list(:,:,ii + fea_group_num*2) = kron(kron(a,bb),cc);
                                witness_list(:,:,ii + fea_group_num*3) = kron(kron(aa,b),cc);
                            end
                        case 12
                            iden = q.si;
                            for ii = 1:fea_group_num
                                a = q.gen_local_random_observe();
                                b = q.gen_local_random_observe();
                                c = q.gen_local_random_observe();
                                aa = q.gen_local_random_observe();
                                bb = q.gen_local_random_observe();
                                cc = q.gen_local_random_observe();

                                witness_list(:,:,ii) = kron(kron(a,b),iden);
                                witness_list(:,:,ii + fea_group_num) = kron(kron(aa,bb),iden);
                                witness_list(:,:,ii + fea_group_num*2) = kron(kron(a,bb),iden);
                                witness_list(:,:,ii + fea_group_num*3) = kron(kron(aa,b),iden);

                                witness_list(:,:,ii + fea_group_num*4) = kron(kron(aa,iden),c);
                                witness_list(:,:,ii + fea_group_num*5) = kron(kron(a,iden),c);
                                witness_list(:,:,ii + fea_group_num*6) = kron(kron(a,iden),cc);
                                witness_list(:,:,ii + fea_group_num*7) = kron(kron(aa,iden),cc); 

                                witness_list(:,:,ii + fea_group_num*8) = kron(kron(iden,b),c);
                                witness_list(:,:,ii + fea_group_num*9) = kron(kron(iden,bb),c);
                                witness_list(:,:,ii + fea_group_num*10) = kron(kron(iden,b),cc);
                                witness_list(:,:,ii + fea_group_num*11) = kron(kron(iden,bb),cc);
                            end
                        case 8
                            for ii = 1:fea_group_num
                                a = q.gen_local_random_observe();
                                b = q.gen_local_random_observe();
                                c = q.gen_local_random_observe();
                                aa = q.gen_local_random_observe();
                                bb = q.gen_local_random_observe();
                                cc = q.gen_local_random_observe();

                                witness_list(:,:,ii) = kron(kron(a,b),c);
                                witness_list(:,:,ii + fea_group_num) = kron(kron(aa,bb),c);
                                witness_list(:,:,ii + fea_group_num*2) = kron(kron(a,bb),c);
                                witness_list(:,:,ii + fea_group_num*3) = kron(kron(aa,b),c);

                                witness_list(:,:,ii) = kron(kron(a,b),cc);
                                witness_list(:,:,ii + fea_group_num) = kron(kron(aa,bb),cc);
                                witness_list(:,:,ii + fea_group_num*2) = kron(kron(a,bb),cc);
                                witness_list(:,:,ii + fea_group_num*3) = kron(kron(aa,b),cc);
                            end
                        case 26
                            for ii = 1:fea_group_num
                                isFirstIden = true;
                                local_wit_1 = q.gen_local_wit(3,isFirstIden);
                                local_wit_2 = q.gen_local_wit(3,isFirstIden);
                                local_wit_3 = q.gen_local_wit(3,isFirstIden);
                                for jj =1:26
                                    p1 = mod(jj,3)+1;
                                    p2 = fix(mod(jj,9)/3)+1;
                                    p3 = fix(jj/9)+1;
                                    witness_list(:,:,ii+(jj-1)*fea_group_num) = ...
                                        kron(kron(local_wit_1(:,:,p1),...
                                        local_wit_2(:,:,p2)),local_wit_3(:,:,p3));             
                                end
                            end  
                        otherwise
                            error('fea_num must be 4, 8, 12 or 26 !!!');
                    end
                 case 4
                     switch fea_num
                         case 80
                             for ii = 1:fea_group_num
                                isFirstIden = true;
                                local_wit_1 = q.gen_local_wit(3,isFirstIden);
                                local_wit_2 = q.gen_local_wit(3,isFirstIden);
                                local_wit_3 = q.gen_local_wit(3,isFirstIden);
                                local_wit_4 = q.gen_local_wit(3,isFirstIden);
                                for jj =1:80
                                    p1 = mod(jj,3)+1;
                                    p2 = fix(mod(jj,9)/3)+1;
                                    p3 = fix(mod(jj,27)/9)+1;
                                    p4 = fix(jj/27)+1;
                                    k1 = kron(local_wit_1(:,:,p1),...
                                        local_wit_2(:,:,p2));
                                    k2 = kron(local_wit_3(:,:,p3),...
                                        local_wit_4(:,:,p4));
                                    witness_list(:,:,ii+(jj-1)*fea_group_num) = ...
                                        kron(k1,k2);            
                                end
                             end
                         case 16
                             for ii = 1:fea_group_num
                                isFirstIden = false;
                                local_wit_1 = q.gen_local_wit(2,isFirstIden);
                                local_wit_2 = q.gen_local_wit(2,isFirstIden);
                                local_wit_3 = q.gen_local_wit(2,isFirstIden);
                                local_wit_4 = q.gen_local_wit(2,isFirstIden);
                                for jj =0:15
                                    p1 = mod(jj,2)+1;
                                    p2 = fix(mod(jj,4)/2)+1;
                                    p3 = fix(mod(jj,8)/4)+1;
                                    p4 = fix(jj/8)+1;
                                    k1 = kron(local_wit_1(:,:,p1),...
                                        local_wit_2(:,:,p2));
                                    k2 = kron(local_wit_3(:,:,p3),...
                                        local_wit_4(:,:,p4));
                                    witness_list(:,:,ii+jj*fea_group_num) = ...
                                        kron(k1,k2);       
                                end
                             end
                         otherwise
                             error('fea_num error')
                     end
                otherwise
                    error('qubit number must be 2 3 or 4!!')
            end
            q.witness_list = witness_list;
            q.data_file = [q.mat_dir , ['CHSH_like_data_for_',num2str(qnum),'_qubits']];
            if fea_group_num>1
                save(q.data_file, 'fea_num');
            end
        end
        
        function witness_list = gen_3_CHSH_like_wit(q)
            disp('generate gen_3_CHSH_like_wit of CHSH-like features');
            
            witness_list = zeros(4,4,9);
            a = q.gen_local_random_observe();
            b = q.gen_local_random_observe();
            aa = q.gen_local_random_observe();
            bb = q.gen_local_random_observe();
            aaa = q.gen_local_random_observe();
            bbb = q.gen_local_random_observe();

            witness_list(:,:,1) = kron(a,b);
            witness_list(:,:,2) = kron(a,bb);
            witness_list(:,:,3) = kron(a,bbb);
            witness_list(:,:,4) = kron(aa,b);
            witness_list(:,:,5) = kron(aa,bb);
            witness_list(:,:,6) = kron(aa,bbb);
            witness_list(:,:,7) = kron(aaa,b);
            witness_list(:,:,8) = kron(aaa,bb);
            witness_list(:,:,9) = kron(aaa,bbb);
            
            q.witness_list = witness_list;
            q.data_file = [q.mat_dir , 'FTI_data_for_2_qubits'];
        end  
        
        function witness_list = gen_multi_wit_bad(q)  % bad
            disp('generate multi qubit features');
            qnum = q.para.qubit_number;
            D = 2^qnum;
            witness_list = nan(D,D,qnum*(qnum-1)/2+1);
            pauli_list = nan(2,2,qnum);
            for i = 1:qnum
                pauli_list(:,:,i) = q.gen_local_random_observe();
            end
            count =1;
            for i = 1:(qnum-1)
                for j = (i+1):qnum
                    flag = zeros(qnum,1);
                    flag(i) = 1;
                    flag(j) = 1;
                    witness = 1;
                    for k = 1:qnum
                        if flag(k)==0
                            witness = kron(witness,eye(2));
                        else
                            witness = kron(witness,pauli_list(:,:,k));
                        end
                    end
                    witness_list(:,:,count) = witness;
                    count = count+1;
                end            
            end
            assert (count == qnum*(qnum-1)/2+1);
            witness = 1;
            for i=1:qnum
                witness = kron(witness,q.gen_local_random_observe());
            end
            witness_list(:,:,end) = witness;
            q.witness_list = witness_list;
           % q.data_file = [q.mat_dir , 'FRO_data_for_multi_qubits'];
        end  
        
        function witness_list = gen_multi_wit(q)  % 为了验证审稿人2的提问
            disp('generate multi qubit features');
            qnum = q.para.qubit_number;
            D = 2^qnum;
            group=2;
            witness_list = nan(D,D,qnum*group);
            for i = 1:qnum*group
                pauli = q.gen_local_random_observe();
                for j = 1:qnum-1
                    pauli = kron(pauli, q.gen_local_random_observe());
                end
                witness_list(:,:,i) = pauli;
            end
            q.witness_list = witness_list;
           % q.data_file = [q.mat_dir , 'FRO_data_for_multi_qubits'];
        end  
        
        function local_wit = gen_local_wit(q,num_of_local_wit,isFirstIden)
            local_wit = zeros(2,2,num_of_local_wit);
            for i = 1:num_of_local_wit
                local_wit(:,:,i) = q.gen_local_random_observe();
            end
            if isFirstIden
                local_wit(:,:,1) = q.si;
            end
        end
        
        function observe = gen_local_random_observe(~)
            r = RandomUnitary(2);
            o1 = r*[0;1];
            o2 = r*[1;0];
            observe = o1*o1'-o2*o2';
        end
        
        function data = extract_data(q,datalen)
            if ~exist('datalen','var')
                datalen = size(q.raw_data);
                datalen = datalen(3);
            end
            q.current_label = q.label(1:datalen);
            disp(['extract_data with length = ',num2str(datalen)]);
            fealen = size(q.witness_list);
            if(length(fealen)>2)
                fealen = fealen(3);
            else
                fealen = 1;
            end
            data = zeros(datalen, fealen);
            % parfor or for
            tic
            parfor d=1:datalen
                for f=1:fealen
                    data(d,f) = trace(q.witness_list(:,:,f)*q.raw_data(:,:,d));
                end
            end
            q.data = real(data);
            toc
        end
        
        function save_data(q,data_file)
            data = q.data;
            label = q.current_label;
            observables = q.witness_list;
            disp('Writing data...');
            if ~exist('data_file','var')
                data_file = q.data_file;
            end
            if exist([data_file,'.mat'],'file')
                save(data_file, 'data','label','observables','-append');
            else
                save(data_file, 'data','label','observables');
            end
        end
        
        function data_name = save_data_for_GHZ_detecting(q)
            data_name = 'Con_Hidden';
            q.data_file = ['mat/',data_name, '_for_', num2str(q.para.qubit_number), '_qubits'];
            q.save_data;       
        end
        
        function save_data_and_p_for_q4(q)
            q.save_data;
            len = q.para.data_length;
            p_random = rand(len,1);
            save(q.data_file,'p_random','-append');
        end
        
        function save_data_and_negativity(q)
            q.save_data;
            negativity = q.negativity;
            save(q.data_file,'negativity','-append');
        end  % for generating gap data
        
        function save_data_for_bisep(q)
            fea_num = q.fea_num;
            disp(['fea_num = ', num2str(fea_num)]);
            q.data_file = ['mat/Bisep_',num2str(q.fea_num)];
            q.save_data;
            save(q.data_file, 'fea_num','-append');
        end
        
        function save_data_for_supp(q,filename,ppt_label,suffix)
            q.data_file = ['mat/',filename];
            if exist('suffix','var')
                q.data_file = [q.data_file,suffix];
            end
            easy_label = q.current_label;
            q.current_label = ppt_label;
            q.save_data;
            save(q.data_file, 'easy_label','-append');
        end
            
    end
       
        

end

