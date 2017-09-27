import colorsys as cs
import random


def hex_to_rgb(value):
    value = value.lstrip('#');
    lv = len(value);
    return tuple(int(value[i:i+lv/3], 16) for i in range(0, lv, lv/3));

def newValue(a, dispA):
    deltaA = random.uniform(-dispA, dispA);
    a = a+deltaA;
    if(a>1):
	a=a-2*abs(deltaA);
    if(a<0):
	a=a+2*abs(deltaA);
    return a;

def newSmartHSVcolor(color, dispH, dispS, dispV):
    r = color[0]/255.0;
    g = color[1]/255.0;
    b = color[2]/255.0;
    
    (h,s,v) = cs.rgb_to_hsv(r, g, b);
#    color = '#%02x%02x%02x' % (color[0], color[1], color[2]);
    
    h=newValue(h,dispH);
    s=newValue(s,dispS);
    v=newValue(v,dispV);

    (r,g,b) = cs.hsv_to_rgb(h, s, v);
    color = '#%02x%02x%02x' % (int(255*r), int(255*g), int(255*b));
    
#    print(h);
#    print(s);
#    print(v);
    return color;

def arrayWithSmartColors(size, dispH, dispS, dispV, color):
    if (color==None):
	r = lambda: random.randint(0,255);
	color = (r(),r(),r());
    else:
	color = hex_to_rgb(color);
    
    array = [];
    for i in range (0, size):
	array.append(newSmartHSVcolor(color, dispH, dispS, dispV));
    return array;
	
	
