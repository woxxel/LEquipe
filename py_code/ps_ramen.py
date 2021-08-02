def ps_ramen(mode):
  scriptSONet([5,1,1,1],[[-1,7],[0],[0],[0]],mode=mode,mode_sim=['pert_spread'],TC=0.2,N=1000)
  scriptSONet([1,6,1,1],[[0],[0,2],[0],[0]],mode=mode,mode_sim=['pert_spread'],TC=0.2,N=1000)
  scriptSONet([1,1,7,1],[[0],[0],[0,6],[0]],mode=mode,mode_sim=['pert_spread'],TC=0.2,N=1000)
  scriptSONet([1,1,1,6],[[0],[0.5],[1],[-0.5,0.5]],mode=mode,mode_sim=['pert_spread'],TC=0.2,N=1000)