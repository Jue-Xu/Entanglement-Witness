classdef GenData < handle
   properties
   end
   
   methods(Static)
       function main_test_run
           % parpool()
           GenData.my_initialize();
           qubit_number = 2;
           raw_data_file = ['mat/raw_data_for_',num2str(qubit_number),'_qubits'];
           data_length = 300000; %300 000 for paper
           min_purity = 0;
           theta_for_val = 10;
           phi_for_val = 60;
           purity_for_val = 100;
           min_theta_for_val = 0;
           para = struct('qubit_number',qubit_number,...
                'raw_data_file',raw_data_file,...
                'data_length',data_length,...
                'min_purity',min_purity,...
                'theta_for_val',theta_for_val,...
                'phi_for_val',phi_for_val,...
                'purity_for_val',purity_for_val,...
                'min_theta_for_val',min_theta_for_val);
        if 1
           delete 'mat\*'
           delete 'npz\*'
           delete 'h5_json\*'
           
           my_rawdata = RawData(para);
           my_rawdata.gen_by_Phi;
        
           my_data = QData(para);
           my_data.parse_raw_data;
           
%          my_data.gen_Tomography_wit;
%          my_data.extract_data();
%          my_data.save_data; 
        
           fea_num = 4; 
           
           my_data.gen_CHSH_like_wit(1,fea_num); % 1 group = no feature selection
       %   my_data.gen_CHSH_wit(); 
           % cancel the comment if you want to "linearly optimize" standard CHSH operators
       
       %   my_data.gen_ShangHaiJiaoTong_wit();   
           % observables for experiment
           my_data.extract_data();
           my_data.save_data_for_GHZ_detecting;
        end
           disp('Generate validation data...')
           my_val_data = ValQCirData(para);
           my_val_data.gen_val_raw_data();
           
           data_name = 'Con_Hidden';
           my_val_data.load_wit_list(data_name);
           my_val_data.extract_data;
           my_val_data.extract_data_by_CHSH;
           disp(['sum of lb_bell = ',num2str(sum(my_val_data.label_bell))]);
           disp(['sum of lb = ',num2str(sum(my_val_data.label))]);
        %   disp(min(my_val_data.label_bell))
        %   disp(max(my_val_data.label_bell))
           my_val_data.save_val_data(strcat('val_CHSH_like_data_for_',num2str(qubit_number),'_qubits'));
           
%            data_name = 'FTI_data';
%            my_val_data.load_wit_list(data_name);
%            my_val_data.extract_data;
%            my_val_data.save_val_data(strcat('val_',data_name));          
       end
       % 2 qubit system where ANN performs better than Bell inequality 
       % (with same measurement resouces)
       
       function q3_BiSep_main
           GenData.my_initialize();
           qubit_number = 3;
           raw_data_file = ['mat/raw_data_for_',...
               num2str(qubit_number),'_qubits'];
           data_length = 200000;  % data_length per group for training  200000 for paper
           para = struct('qubit_number',qubit_number,...
                'data_length',data_length,...
                'raw_data_file',raw_data_file);           
           if 1
               delete 'mat\*'
               delete 'npz\*'
               delete 'h5_json\*'
               my_rawdata = RawData(para);        
           %    my_rawdata.gen_by_GHZ_W_Channel;
           
               my_rawdata.gen_by_BiSep;
           end
           my_data = QData(para);
           my_data.parse_raw_data;
           
           my_data.gen_Tomography_wit;
           my_data.extract_data;   
           my_data.save_data('mat\Bisep_63.mat');
           
           fea_num = 26;
           my_data.gen_CHSH_like_wit(1,fea_num);  % 1 group is enough
           my_data.extract_data;   
           my_data.save_data('mat\Bisep_26.mat');
           
           fea_num = 12;
           my_data.gen_CHSH_like_wit(1,fea_num);  % 1 group is enough
           my_data.extract_data;   
           my_data.save_data('mat\Bisep_12.mat');
           
           fea_num = 8;
           my_data.gen_CHSH_like_wit(1,fea_num);  % 1 group is enough
           my_data.extract_data;   
           my_data.save_data('mat\Bisep_8.mat');
           
           fea_num = 4;
           my_data.gen_CHSH_like_wit(1,fea_num);  % 1 group is enough
           my_data.extract_data;   
           my_data.save_data('mat\Bisep_4.mat'); 
       end
       
       function q4_main_paper   % paper for 4 qubit
           GenData.my_initialize();
           qubit_number = 4;
           data_length = 170000;  % 170000 data_length per group for training BEll-like
           raw_data_file = 'mat/raw_data_for_4_qubits';
           div = 7;
           max_noise = 1;
           para = struct('qubit_number',qubit_number,...
                'data_length',data_length,...
                'raw_data_file',raw_data_file,...
                'max_noise',max_noise,...
                'div',div);  
           if 1
               delete 'mat\*'
               delete 'npz\*'
               delete 'h5_json\*'
               my_rawdata = RawData(para);
               
          %     raw_data = my_rawdata.gen_by_EntanglementWithNoise;
          %     lb = RawData.check_by_PPT(raw_data);
          %     sum(lb)/length(lb)
          %     my_rawdata.gen_by_Division;
               
               my_rawdata.gen_by_PureState;  %pure sep sep ---- The order must be Right!!
               my_rawdata.gen_by_Division;         
               my_rawdata.gen_by_Division;
           end
           my_data = QData(para);
           my_data.parse_raw_data;
            
           fea_num = 80;
           my_data.gen_CHSH_like_wit(1,fea_num);
           my_data.extract_data;         
           my_data.save_data_and_p_for_q4;
       end
      
       function q4_main   % for 4 qubit (just test, not in paper)
           GenData.my_initialize();
           qubit_number = 4;
           data_length = 70000;
           max_noise = 1;
           raw_data_file = 'mat/raw_data_for_4_qubits';
           div = 7;
           para = struct('qubit_number',qubit_number,...
                'data_length',data_length,...
                'raw_data_file',raw_data_file,...
                'max_noise',max_noise,...
                'div',div);
           if 1
               delete 'mat\*'
               delete 'npz\*'
               delete 'h5_json\*'
               my_rawdata = RawData(para);
               
               my_rawdata.gen_by_EntanglementWithNoise;
               my_rawdata.gen_by_Division;              
           
               my_data = QData(para);
               my_data.parse_raw_data;
               lb = RawData.check_by_PPT(my_data.raw_data); 
               fea_num = 80;
               observables = my_data.gen_CHSH_like_wit(1,fea_num);
              % my_data.gen_Tomography_wit();
               my_data.extract_data;
               my_data.save_data_for_supp('q4_80',lb);
            end
               para.data_length = 0.4*para.data_length;  
               para.raw_data_file = [para.raw_data_file,'_test'];
               my_rawdata = RawData(para);
               my_rawdata.gen_by_EntanglementWithNoise;
               my_rawdata.gen_by_Division;           
           
               my_data = QData(para);
               my_data.parse_raw_data;
               lb = RawData.check_by_PPT(my_data.raw_data);
               %%%%% 教训：check_by_PPT的错误，min() 需要改成 min(real())，已改正 %%%
               %%% 话说Matlab 2016和 2017 在min(-3+0.0i,2+0i)的结果似乎不一样，未来需验证

               
               my_data.witness_list = observables;
           %    q.data_file = [q.mat_dir , ['CHSH_like_data_for_4_qubits']];
               my_data.extract_data; 
               my_data.save_data_for_supp('q4_80',lb,'_test');
       end
       
       function main_general_detection   
           GenData.my_initialize();
           raw_data_file = 'mat/raw_data_for_gap';
           qubit_number = 2;
           data_length = 300000;  %3,000,000 for train, 300,000 for test
           para = struct('qubit_number',qubit_number,...
                'raw_data_file',raw_data_file,...
                'data_length',data_length,...
                'fill',false); 
        istrain = 0;
        if 0
           delete 'mat\*'
           delete 'npz\*'
           delete 'h5_json\*'
        end
        if ~istrain
           para.raw_data_file = 'mat/RawGapData_test';
        end
           my_rawdata = GapRawData(para);
           my_rawdata.gen_by_Gap;
        
           my_data = QData(para);
           my_data.parse_raw_data(1); %1 = load negativity

           my_data.gen_Tomography_wit;
           my_data.extract_data();
        if ~istrain
           my_data.data_file = 'mat/Gap_data_for_test.mat';
        end
           my_data.save_data();
           negativity = my_data.negativity;
           save(my_data.data_file,'negativity','-append');
       end
       % 2 qubit system for general entanglement detection

       function main_UPB   % generate bound entangled states
           GenData.my_initialize();
           qubit_number = 3;
           data_length = 200000; % 200000 for paper
           raw_data_file = 'mat/raw_data_for_3_qubits';
           div = 17;
           para = struct('qubit_number',qubit_number,...
                'data_length',data_length,...
                'raw_data_file',raw_data_file,...;
                'div',17);
           if 0
               delete 'mat\*'
               delete 'npz\*'
               delete 'h5_json\*'
               my_rawdata = RawData(para);
               
               my_rawdata.gen_Bisep_for_UPB;
               my_rawdata.gen_by_UPB;               
           end
           my_data = QData(para);
           my_data.parse_raw_data;

           my_data.gen_CHSH_like_wit(1,4);
           my_data.extract_data;      
           my_data.save_data('mat\UPB_4.mat');
            
           my_data.gen_CHSH_like_wit(1,8);
           my_data.extract_data;      
           my_data.save_data('mat\UPB_8.mat');
           
           my_data.gen_CHSH_like_wit(1,12);
           my_data.extract_data;      
           my_data.save_data('mat\UPB_12.mat');
           
           my_data.gen_CHSH_like_wit(1,26);
           my_data.extract_data;      
           my_data.save_data('mat\UPB_26.mat');
           
           my_data.gen_Tomography_wit;
           my_data.extract_data;      
           my_data.save_data('mat\UPB_63.mat');
       end
       
       function my_initialize()
           p1 = mfilename('fullpath');
           i=strfind(p1,'\');
           p1=p1(1:i(end));
           cd(p1);
           rng(0);   
       end
       
       function main_diri()
           GenData.my_initialize();
           raw_data_file = 'mat/raw_data_for_diri';
           qubit_number = 3;
           data_length = 10000;  %700000!!
           para = struct('qubit_number',qubit_number,...
                'raw_data_file',raw_data_file,...
                'data_length',data_length,...
                'fill',false,...
                'min_neg',-0.01,...
                'max_neg',0.01);         
        if 1
           delete 'mat\*'
           delete 'npz\*'
           delete 'h5_json\*'
           
           my_rawdata = GapRawData(para);
           my_rawdata.gen_by_Diri;
        end
           my_data = QData(para);
           my_data.parse_raw_data(1); %1 means load negativity

           my_data.gen_Tomography_wit;
           my_data.extract_data();
           my_data.save_data('mat\FTI_data_for_2_qubits.mat')
       end
       % generate random states (diri distribution) in order to 
       % compare with ZengBei's work 
           
       function q_multi_for_exp   % multi-qubit system for test scaling cases
           GenData.my_initialize();
           data_length = 500;
           div = 7;
           para = struct('data_length',data_length,...
                'div',div);
           if 1
               delete 'mat\*'
               delete 'npz\*'
               delete 'h5_json\*'
           end
                           
           for qubit_number = 2:8
               for repeat = 1:10
                   data_name = [num2str(qubit_number), '_qubits_',num2str(repeat),'repeat'];
                   raw_data_file = ['mat/raw_data_for_',data_name];
                   para.raw_data_file = raw_data_file;
                   para.qubit_number=qubit_number;

                   my_rawdata = RawData(para);
                   my_rawdata.gen_GHZ;
                   my_rawdata.gen_by_Division;
                   my_data = QData(para);
                   my_data.parse_raw_data;
                   my_data.gen_multi_wit();
                   my_data.extract_data();
                   my_data.save_data(['mat/data_for_',data_name]);
               end
           end
       end
   end
end

