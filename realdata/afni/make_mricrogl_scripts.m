FDRmethod = 'bky2006';
switch lower(FDRmethod)
    case {'bh1995','bh'}
        fdrstr = 'bh';
    case {'bky2006','bky'}
        fdrstr = 'bky';
end

txt = ['import gl\n' ...
'gl.resetdefaults()\n' ...
'gl.bmpzoom(2)\n' ...
'gl.loadimage("registered.nii.gz")\n' ...
'gl.overlayload("stats.FT.masked.nii.gz")\n' ...
'gl.minmax(1, %g, 15)\n' ...
'gl.colorfromzero(1,1)\n' ...
'gl.opacity(1,70)\n' ...
'gl.colorname (1,"6warm")\n' ...
'gl.overlayload("stats.FT.masked.nii.gz")\n' ...
'gl.minmax(2, -15, %g)\n' ...
'gl.colorfromzero(1,1)\n' ...
'gl.opacity(1,70)\n' ...
'gl.colorname (2,"5winter")\n' ...
'gl.overlaymaskwithbackground(0)\n' ...
'gl.mosaic("A L- V -0.1 -32 -16 0 8 16 32 40 50 S X R 0")\n' ...
'gl.savebmp("./%s_%s.png")\n' ...
'exit()\n'];

Jfiles = dir('*.json');
for j = 1:numel(Jfiles)
    [~,jname,~] = fileparts(Jfiles(j).name);
    method = split(jname,'_'); method = method{end};
    J = readjson(Jfiles(j).name);
    if isempty(J.zposthr), J.zposthr = 7; end
    if isempty(J.znegthr), J.znegthr = 7; end
    if strcmpi(J.zposthr,'does not apply'), J.zposthr = 7; end
    if strcmpi(J.znegthr,'does not apply'), J.znegthr = 7; end
    fid = fopen(sprintf('make_fig_%s_%s.py',fdrstr,method),'w');
    fprintf(fid,txt,J.zposthr,-J.znegthr,fdrstr,method);
    fclose(fid);
end