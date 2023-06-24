function J = readjson(jsonfile)
fid  = fopen(jsonfile);
Jstr = fread(fid,inf);
Jstr = char(Jstr');
fclose(fid);
J = jsondecode(Jstr);