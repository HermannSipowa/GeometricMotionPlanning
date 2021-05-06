Gmailserver = 'smtp.gmail.com';
Myaddress = 'kaptuihremann@gmail.com';
Mypassword = 'hil1de9gar9de4';
Emailsubject = 'Matlab Code on Fortuna';
Emailmsg = 'Hey man, Your code is done running. See the attached files below. Best, Fortuna.';
Attachment0 = 'MatlabOutput.txt';
Attachment1 = 'q_t0.mat';
Attachment2 = 'q_tf.mat';
Attachment3 = 'sol.mat';
Attachment4 = 'Integrated_Trajectory_6DoF2.png';
Attachment5 = 'Quaternion_Tracking_6DoF.png';


setpref('Internet','E_mail',Myaddress);
setpref('Internet','SMTP_Server',Gmailserver);
setpref('Internet','SMTP_Username',Myaddress);
setpref('Internet','SMTP_Password',Mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
sendmail(Myaddress,Emailsubject,Emailmsg,{Attachment0,Attachment1,...
    Attachment2,Attachment3});