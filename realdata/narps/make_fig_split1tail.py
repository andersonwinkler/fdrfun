import gl
'''
import json

# ============================================================================
def readjson(jsonfile):
    """Read a json file to dict."""
    with open(jsonfile, 'r') as fp:
        J = json.load(fp)
    return J

def writejson(J, jsonfile, epoch=None): 
    """Write a dict to json file."""
    with open(jsonfile, 'w') as fp:
        json.dump(J, fp, indent=2)
    if epoch is not None:
        os.utime(jsonfile, (epoch,epoch))
    return
# ============================================================================
'''
gl.resetdefaults()
gl.bmpzoom(2)
gl.loadimage('MNI152_T1_2mm_brain.nii')
#open overlay: show positive regions
gl.overlayload('narps-4965_9U7M-hypo1_split1tail.nii')
gl.minmax(1, 2, 7)
gl.opacity(1,70)
gl.colorname (1,"6warm")
#open overlay: show negative regions
gl.overlayload('narps-4965_9U7M-hypo1_split1tail.nii')
gl.minmax(2, -2, -7)
gl.opacity(1,70)
gl.colorname (2,"5winter")
gl.overlaymaskwithbackground(0)
#to mask with background:
# gl.overlaymaskwithbackground(1)
gl.mosaic("A L- V -0.1 -32 -16 0 8 16 32 40 50 S X R 0")
gl.savebmp("./narps.png")
