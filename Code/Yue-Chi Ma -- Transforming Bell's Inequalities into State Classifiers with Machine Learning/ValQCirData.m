classdef ValQCirData < handle
    %UNTITLED 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        para
        data
        raw_data
        label
        label_bell
        witness_list
        theta
        phi_num
        purity
        
        si = eye(2);
        sx=[0,1;1,0];
        sy=[0,-1i;1i,0];
        sz=[1,0;0,-1];
    end
    
    methods
        function q = ValQCirData(para)
            q.para = para;
        end
        
        function load_wit_list(q,data_name)
            wl = load(['mat/' , data_name,'_for_',num2str(q.para.qubit_number),'_qubits']);
            wl = wl.observables;
            q.witness_list=wl;
        end
        
        function data = extract_data(q)
            disp('extracting data...');
            datalen = size(q.raw_data);
            datalen = datalen(3);
            fealen = size(q.witness_list);
            if(length(fealen)>2) %check if witness_list only has one operator(feature)
                fealen = fealen(3);
            else
                fealen = 1;
            end
            data = zeros(datalen, fealen);
            disp(['validation datalen = ',num2str(datalen)]);
            
            
            parfor d=1:datalen
                for f = 1:fealen
                    data(d,f) = trace(q.witness_list(:,:,f)*q.raw_data(:,:,d));
                end
            end
            q.data = real(data);
        end
        
        function label = extract_data_by_CHSH(q)
            datalen = size(q.raw_data);
            datalen = datalen(3);
            qdata = QData(q.para);
            witness_list = qdata.gen_CHSH_wit();
            fealen = 4;
            assert(length(witness_list) == fealen)
            data = zeros(datalen, fealen);
            label = zeros(datalen,1);
            disp(['validation datalen = ',num2str(datalen)]);
            parfor d=1:datalen
                label(d) = ...
                    trace(q.raw_data(:,:,d)*witness_list(:,:,1)-...
                    q.raw_data(:,:,d)*witness_list(:,:,2)+...
                    q.raw_data(:,:,d)*witness_list(:,:,3)+...
                    q.raw_data(:,:,d)*witness_list(:,:,4));
            end
            q.label_bell = abs(label)<=2;
        end
        
        function raw_data = gen_val_raw_data(q)
            p = q.para;
            qubit_number = p.qubit_number;
            D = 2^qubit_number;
            theta_num = p.theta_for_val;
            phi_num = p.phi_for_val;
            purity_num = p.purity_for_val;
            min_theta = p.min_theta_for_val;
            theta = linspace(min_theta,pi/4,theta_num);
            phi = linspace(0,2*pi,phi_num+1);
            phi = phi(1:end-1);
            cos_theta = cos(theta);
            sin_theta = sin(theta);
            exp_i_phi = exp(1i*phi);
            purity = linspace(0,1,purity_num);
            each_tp_num = theta_num*phi_num;
            max_mixure = repmat(eye(D)/D,1,1,each_tp_num);
            tmp_raw_data = zeros(D,D,each_tp_num);
            lb = zeros(each_tp_num*purity_num,1);
            n=1;
%             switch qubit_number
%                 case 2
%                     s_state = [1,0,0,1]';
%                 otherwise
%                     error('qubit_number must be 2 !!');
%             end
            
            for t = 1:theta_num
                for my_phi = 1:phi_num
                    state = [0,cos_theta(t),exp_i_phi(my_phi)*sin_theta(t),0]';
                    tmp_raw_data(:,:,n) = state * state';
                    n = n+1;
                end
            end
            
            raw_data = zeros(D,D,each_tp_num*purity_num);

            for pu = 1:purity_num
                raw_data(:,:,(pu-1)*each_tp_num+1:pu*each_tp_num)...
                    = purity(pu)*tmp_raw_data + (1-purity(pu))*max_mixure;
            end
            
            % sys = ones(1,qubit_number)*2;
            disp('generate val_data label...');
            
            parfor l=1:each_tp_num*purity_num
                lb(l) = IsPPT(raw_data(:,:,l));
                
             %   lb(l) = min(real(eig(...
              %              PartialTranspose(raw_data(:,:,l),1,sys))))>=0
            end
            
            q.raw_data = raw_data;
            q.label = lb;
            q.purity = purity;
            q.theta = theta;
            q.phi_num = phi_num;
        end
        
        function raw_data = gen_val_raw_data_for_3ran_U(q)
            % for GHZ-class state validation
            
            D = 8;
            theta_num = q.para.theta_for_val;
            q.phi_num = q.para.phi_for_val;
            U_list = zeros(D,D,q.phi_num);
            for ph = 1:q.phi_num
                U_list(:,:,ph) = kron(kron(RandomUnitary(2),...
                    RandomUnitary(2)),RandomUnitary(2));
            end 
            purity_num = q.para.purity_for_val;
            q.theta = linspace(0.05,pi/4,theta_num);

            cos_theta = cos(q.theta);
            sin_theta = sin(q.theta);
            
            q.purity = linspace(0,1,purity_num);
            each_tp_num = theta_num*q.phi_num;
            max_mixure = repmat(eye(D)/D,1,1,each_tp_num);
            tmp_raw_data = zeros(D,D,each_tp_num);
            lb = zeros(each_tp_num*purity_num,1);
        %    trace_rem = zeros(each_tp_num,1);
            ori_state = [1,0,0,0,0,0,0,1]'/sqrt(2);
            for t = 1:theta_num
                state = ori_state;
                state(1) = cos_theta(t);
                state(end) = sin_theta(t);
                for ph = 1:q.phi_num
                    U = kron(kron(RandomUnitary(2),...
                        RandomUnitary(2)),RandomUnitary(2));
                    % new_state = U_list(:,:,ph) * state;
                    new_state = U * state;
                    tmp_raw_data(:,:,ph + (t-1)*q.phi_num) = new_state * new_state';
                end
            end
            raw_data = zeros(D,D,each_tp_num*purity_num);

            for pu = 1:purity_num
                raw_data(:,:,(pu-1)*each_tp_num+1:pu*each_tp_num)...
                    = q.purity(pu)*tmp_raw_data + (1-q.purity(pu))*max_mixure;
             %   lb((pu-1)*each_tp_num+1:pu*each_tp_num)...
             %       = q.purity(pu)*trace_rem + (1-q.purity(pu))*max_mix_trace;
            end
       %     lb = lb>0;
            
            parfor l=1:each_tp_num*purity_num
                lb(l) = min(eig(...
                    PartialTranspose(raw_data(:,:,l),1,[2,2,2])))>=0;
            end
            q.raw_data = raw_data;
            q.label = lb;
        end
        
        function save_val_data(q, data_name)
            disp(['Writing ',data_name]);
            data = q.data;
            label = q.label;
            label_bell = q.label_bell;
            theta_list = q.theta;
            purity_list = q.purity;
            phi_num = q.phi_num;
            save(strcat('mat/',data_name),'data','label',...
                'theta_list','purity_list','label_bell','phi_num');
        end   
    end
end

