% dianpd.mp -- do not edit, generated automatically by presentation.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 80*1.0pt#, 50*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
straight_galleries;
vleft (__i1, __i2, __i3, __i4, __i5, __i6);
vright (__o1, __o2, __o3, __o4, __o5, __o6);
vconnect ("plain", __i1, __v1, __o1);
vconnect ("plain", __i2, __v2, __o2);
vconnect ("plain", __i3, __v3, __o3);
vconnect ("plain,tension=0.4", __i4, __v4);
vconnect ("plain,tension=0.4", __v5, __o4);
vconnect ("fermion,label=$\vec {k}$,label.side=left", __v4, __v6);
vconnect ("fermion,label=$\vec {k}$, label.side=left", __v7, __v5);
vdot (__v4, __v5);
vforce (((0.3w,1.0h)), __v6);
vforce (((0.7w,1.0h)), __v7);
vforce (((0.3w,0.6h)), __v4);
vforce (((0.7w,0.6h)), __v5);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
