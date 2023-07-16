import gl
gl.resetdefaults()
gl.bmpzoom(2)
gl.loadimage("registered.nii")
gl.overlayload("stats.FT.masked.nii")
gl.minmax(1, 2.31328, 15)
gl.opacity(1,70)
gl.colorname (1,"6warm")
gl.overlayload("stats.FT.masked.nii")
gl.minmax(2, -15, -7)
gl.opacity(1,70)
gl.colorname (2,"5winter")
gl.overlaymaskwithbackground(0)
gl.mosaic("A L- V -0.1 -32 -16 0 8 16 32 40 50 S X R 0")
gl.savebmp("./bky_canonical.png")
exit()
