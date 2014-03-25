This code implements the local community detection methods used in 
    Jeub, L. G. S., Balachandran, P., Porter, M. A., Mucha, P. J., 
    & Mahoney, M. W. (2014, March 15). 
    Think Locally, Act Locally: The Detection of Small, Medium-Sized, 
    and Large Communities in Large Networks. arXiv:1403.3795 [cs.SI]

The main interface is the NCP function which computes local and global 
NCPs using the AclCut, MovCut, or EgoNet methods described in the paper. 

To create a conductance ratio profile (CRP), first compute the network 
community profile using the NCP function, storing both the conductance 
values and communities. Then compute the internal conductance for 
these communities using the InternalConductance function. One then 
obtains the conductance ratio profile by dividing the NCP by the 
internal conductance, e.g.

[conductance,S]=NCP(A,'ACL');
[internal_conductance]=InternalConductance(A,S);

CRP=conductance./internal_conductance

This code is licensed under a FreeBSD license (see License.txt), except
for lobpcg.m by A.V. Knyazev, which is provided under the GNU LGPL ver 2.1.
