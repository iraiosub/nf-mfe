"""Map observed RNAduplex windows before counting GC at paired positions."""
import numpy as np

_BASE = np.zeros(256,dtype=np.int8)
for code,base in enumerate('ACGU',1):
    _BASE[ord(base)]=code
_PAIR = np.zeros((5,5),dtype=bool)
for l,r in [('A','U'),('U','A'),('C','G'),('G','C'),('G','U'),('U','G')]:
    _PAIR[_BASE[ord(l)],_BASE[ord(r)]]=True


def gc_at_window(left,right,lmask,rmask,lstart,rstart):
    """Both partners contribute one nucleotide; unpaired dots never contribute."""
    if lstart<0 or rstart<0 or lstart+len(lmask)>len(left) or rstart+len(rmask)>len(right):
        raise ValueError('Duplex window outside arm sequence')
    paired=''.join(base for base,mark in zip(left[lstart:lstart+len(lmask)],lmask) if mark=='(')
    paired+=''.join(base for base,mark in zip(right[rstart:rstart+len(rmask)],rmask) if mark==')')
    if not paired:
        return np.nan, 'zero_pairs'
    if any(base not in 'ACGU' for base in paired):
        return np.nan, 'ambiguous_paired_nucleotide'
    return 100.0*sum(base in 'GC' for base in paired)/len(paired),'ok'


def paired_gc_percent(task, resolver=None):
    """Return GC%, mapping status, and actual paired-nucleotide count.

    A unique canonical alignment fixes the window exactly. If all possible
    canonical alignments give the same GC count, GC is likewise exact even
    when the offsets differ. Otherwise RNAduplex recovers i/j, and its saved
    structure and energy must agree before its coordinates are used.
    """
    left,right,structure,mfe=task
    left=str(left).strip().upper().replace('T','U')
    right=str(right).strip().upper().replace('T','U')
    lmask,rmask=str(structure).split('&')
    lpos=np.array([i for i,c in enumerate(lmask) if c=='('],dtype=int)
    rpos=np.array([i for i,c in enumerate(rmask) if c==')'][::-1],dtype=int)
    if len(lpos)!=len(rpos):
        raise ValueError('Unbalanced paired positions')
    count=2*len(lpos)
    if count==0:
        return np.nan,'zero_pairs',0
    nl,nr=len(left)-len(lmask)+1,len(right)-len(rmask)+1
    if nl<=0 or nr<=0:
        raise ValueError('Dot-bracket window longer than sequence')
    lb=_BASE[np.frombuffer(left.encode('ascii'),dtype=np.uint8)]
    rb=_BASE[np.frombuffer(right.encode('ascii'),dtype=np.uint8)]
    ls=np.repeat(np.arange(nl),nr);rs=np.tile(np.arange(nr),nl)
    for lp,rp in zip(lpos,rpos):
        compatible=_PAIR[lb[ls+lp],rb[rs+rp]]
        ls,rs=ls[compatible],rs[compatible]
        if not len(ls):
            break
        if len(ls)==1:
            if _PAIR[lb[ls[0]+lpos],rb[rs[0]+rpos]].all():
                value,status=gc_at_window(left,right,lmask,rmask,int(ls[0]),int(rs[0]))
                return value,'unique_alignment' if status=='ok' else status,count
            ls=ls[:0];rs=rs[:0]
            break
    if len(ls)>1:
        gc=((lb[ls[:,None]+lpos]==2)|(lb[ls[:,None]+lpos]==3)).sum(axis=1)
        gc+=((rb[rs[:,None]+rpos]==2)|(rb[rs[:,None]+rpos]==3)).sum(axis=1)
        if (gc==gc[0]).all():
            return float(100.0*gc[0]/count),'gc_invariant_alignments',count
    if resolver is None:
        import RNA
        resolver=RNA.duplexfold
    duplex=resolver(left,right)
    if duplex.structure!=structure or not np.isclose(duplex.energy,mfe,rtol=0,atol=1e-4):
        raise ValueError('Cannot recover saved duplex: structure/energy mismatch')
    value,status=gc_at_window(left,right,lmask,rmask,int(duplex.i)-len(lmask),int(duplex.j)-1)
    return value,'verified_duplex_coordinates' if status=='ok' else status,count
