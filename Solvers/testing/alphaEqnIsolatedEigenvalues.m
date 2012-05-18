LL=U;
LL.internal=U.internal+alphag.internal.*(-0.282.*(1-alphag.internal./(alphag.internal+1000.*(1-alphag.internal)))-(0.282.*(1-alphag.internal))./(alphag.internal+1000.*(1-alphag.internal)))+0.282.*(1-alphag.internal).*(1-alphag.internal./(alphag.internal+1000.*(1-alphag.internal)));
LL.left.type='V';
LL.left.value=0;
LL.right.type='V';
LL.right.value=0;
LL=setBC(LL,constField(0,N),xC,xF,g);