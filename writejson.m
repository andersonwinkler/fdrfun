function writejson(J,jsonfile)
Jstr = jsonencode(J,'PrettyPrint',true);
fid  = fopen(jsonfile,'w');
fprintf(fid,Jstr);
fclose(fid);
