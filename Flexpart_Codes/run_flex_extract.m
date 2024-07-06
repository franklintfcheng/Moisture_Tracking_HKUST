%%

error('This code is aborted. Use the python script "python_run_flex_extract.py" instead -> More stable.')

% %======= Assign the work to your ECMWF farmers! =======%
% flex_extract_path = '/disk/r059/tfchengac/flexpart_v10.4/preprocess/flex_extract_v7.1.2/Run/';
% slct_year_list = 2017;  st_month = 2; st_day = 8; ed_month = 3;
% 
% % flex_extract_path = '/disk/r034/tfchengac/flex_extract_farmers/farmer_1/Run/';
% % slct_year_list = 2018;  st_month = 1; ed_month = 3;
% 
% % flex_extract_path = '/disk/r034/tfchengac/flex_extract_farmers/farmer_2/Run/';
% % slct_year_list = 2017;  st_month = 10; ed_month = 12;
% 
% % flex_extract_path = '/disk/r034/tfchengac/flex_extract_farmers/farmer_3/Run/';
% % slct_year_list = 2018; st_month = 10; ed_month = 12;
% % 
% % farmer_id = 1; %/ assign farmers 101-112 to download data for each month
% % farmer_id_str = sprintf('%02d',farmer_id); %/ convert int into str with leading zeros
% % 
% % flex_extract_path = strcat('/disk/r034/ldaiad/flex_extract_farmers/farmer_1',farmer_id_str,'/Run/');
% % slct_year_list = 2016:-1:2010;  st_month = farmer_id; ed_month = farmer_id;
% 
% dates = [];
% for slct_year = slct_year_list
% 
%     date_st_dt = datetime(slct_year, st_month,  st_day, 'format','yyyyMMdd');
%     date_ed_dt = datetime(slct_year, ed_month,  eomday(slct_year,ed_month), 'format','yyyyMMdd'); %/ use eomday() to find the last day of the month
% 
% %     date_st_dt = datetime(slct_year, st_month,  1, 'format','yyyyMMdd');
% %     date_ed_dt = datetime(slct_year, ed_month,  eomday(slct_year,ed_month), 'format','yyyyMMdd'); %/ use eomday() to find the last day of the month
% 
%     dates_dt = [date_st_dt:date_ed_dt]';
%     dates = cat(1, dates, yyyymmdd(dates_dt));
% end
% 
% 
% for t = 1:length(dates)
%     tic
%     disp(flex_extract_path)
%     disp(dates(t))
%     status = [];
%     cd([flex_extract_path, 'Control']);
%     
%     CONTROL_filename = 'CONTROL_EA5.global.matlab';
%     
%     CONTROL_content = {
%                     ['START_DATE ', num2str(dates(t))],...
%                     'DTIME 1',...
%                     'TYPE AN AN AN AN AN AN AN AN',...
%                     'TIME 00 03 06 09 12 15 18 21',...
%                     'STEP 00 00 00 00 00 00 00 00',...
%                     'ACCTYPE FC',...
%                     'ACCTIME 06/18',...
%                     'ACCMAXSTEP 12',...
%                     'CLASS EA',...
%                     'STREAM OPER',...
%                     'DATASET ERA5',...
%                     'GRID 1',...
%                     'LEFT -179.',...
%                     'LOWER -90.',...
%                     'UPPER 90.',...
%                     'RIGHT 180.',...
%                     'LEVELIST 1/to/137',...
%                     'RESOL 159',...
%                     'ETA 1',...
%                     'CWC 1',...
%                     'RRINT 1',...
%                     'FORMAT GRIB2',...
%                     'PREFIX EA',...
%                     'ECTRANS 0'};
% 
%     fid=fopen(CONTROL_filename,'w');
%     for l = 1:length(CONTROL_content)
%         fprintf(fid, strcat(CONTROL_content{l}, '\n'));
%     end
%     fclose(fid);
%     
%     cd(flex_extract_path);
% %     system('ls');
%     logfilename = ['EA', num2str(dates(t)), '.log'];
% 
%     command = ['nohup ./run.sh > ', logfilename];
%     disp(command)
%     [status, cmdout] = system(command);
%     disp(cmdout)
%     disp(status)
%     toc
% end

