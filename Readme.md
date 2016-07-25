# LocalCommunities

This code implements the local community detection methods used in
 
1.	Jeub, L. G. S., Balachandran, P., Porter, M. A., Mucha, P. J., & Mahoney, M. W. (2015). 
	Think locally, act locally: Detection of small, medium-sized, and large communities in large networks. 
	Physical Review E, 91(1), 012821. 
	http://doi.org/10.1103/PhysRevE.91.012821

2. 	Jeub, L. G. S., Mahoney, M. W., Mucha, P. J., & Porter, M. A. (2015). 
	A Local Perspective on Community Structure in Multilayer Networks. 
	arXiv:1510.05185 [cs.SI].
	http://arxiv.org/abs/1510.05185

## Usage

The main interface is the `NCP` function which computes local and global NCPs using the `ACLcut`, `MOVcut`, or `EGOcut` methods described in [1]. For multilayer networks `Multilayer/multilayerNCP` provides a convenience wrapper around `NCP` that sets up `NCP` with appropriate options for the ‘classical’ or ‘relaxed’ random walk and allows for sampling using either ‘physical’ or ‘state’ nodes.

To create a conductance ratio profile (CRP), first compute the network community profile using the `NCP` function, storing both the conductance values and communities. Then compute the internal conductance for these communities using the `InternalConductance` function. One then obtains the conductance ratio profile by dividing the NCP by the internal conductance, e.g.

    [conductance,S]=NCP(A,'ACL');
    [internal_conductance]=InternalConductance(A,S);
    
    CRP=conductance./internal_conductance

More detailed documentation of the different functions is available using `help LocalCommunities` or `doc LocalCommunities` in Matlab.

## Installation

To use the code, make sure that the ‘LocalCommunities’ folder and its subfolders are on your MATLAB path. Next, compile the `APPR` mex function by running `compile_APPR_mex` in the ‘ACLcut’ folder.

Some functions require the ‘matlab_bgl’ library which is available from http://www.mathworks.co.uk/matlabcentral/fileexchange/10922-matlabbgl

## License and Citation

This code is licensed under a FreeBSD license (see License.txt), except
for lobpcg.m by A.V. Knyazev, which is provided under the GNU LGPL ver 2.1.

When using the code, please cite the relevant papers and include a pointer to the ‘LocalCommunities’ Github page https://github.com/LJeub/LocalCommunities . If you want to cite the code directly, include a citation similar to:

* Jeub, L. G. S. (2014—2016). LocalCommunities. https://github.com/LJeub/LocalCommunities

