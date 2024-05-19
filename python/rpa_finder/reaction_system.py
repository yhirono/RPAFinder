import os
import numpy as np
import sympy as sp
from scipy.linalg import null_space
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import libsbml
from .utils import extract_species, from_str_to_list, flatten_list, parse_string_to_reaction_list
from .linalg import (
    proj_mat, proj_mat_comp, dimker, dimcoker, 
    get_intersection, get_intersection_symb, 
    get_orthogonal_complement, get_orthogonal_complement_symb
)


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})


def _sort_rows_and_columns_by_true_count(boolean_matrix):
    # Count the number of True values in each row and column
    row_true_counts = np.sum(boolean_matrix, axis=1)
    col_true_counts = np.sum(boolean_matrix, axis=0)

    # Sort the rows based on the number of True values in decreasing order
    sorted_row_indices = np.argsort(row_true_counts)[::-1]

    # Sort the columns based on the number of True values in decreasing order
    sorted_col_indices = np.argsort(col_true_counts)[::-1]

    # Reorder the rows and columns of the boolean matrix
    sorted_matrix = boolean_matrix[sorted_row_indices][:, sorted_col_indices]

    return sorted_matrix, sorted_row_indices, sorted_col_indices


def _contains_all(list1, list2):
    return all(element in list1 for element in list2)


def _randomize_nonzero_entries(matrix):
    nonzero_mask = matrix != 0
    random_numbers = np.random.uniform(1., 1.5, size=matrix.shape)
    return np.where(nonzero_mask, random_numbers, 0)


def _make_hashable(obj):
    if isinstance(obj, (tuple, list)):
        return tuple((_make_hashable(e) for e in obj))
    elif isinstance(obj, np.ndarray):
        return tuple(obj.tolist())
    else:
        return obj


def _remove_duplicates(lst):
    seen = set()
    new_lst = []
    for item in lst:
        hashable_item = _make_hashable(item)
        if hashable_item not in seen:
            seen.add(hashable_item)
            new_lst.append(item)
    return new_lst


def _get_indices_of_element(the_list, element):
    return [index for index, item in enumerate(the_list) if item == element]


def from_lbs_to_bs(lbs):
    """
    Convert a labeled buffering structure to a buffering structure (i.e., a subnetwork specified by species and reactions).

    Args:
        lbs (list): A labeled buffering structure represented by a list of four lists.

    Returns:
        list: A list of lists of species and reactions, representing a subnetwork.
    """
    return [lbs[1], sorted(lbs[2] + lbs[3])]


def _nonzero_with_error(arr, ep=1e-5):
    return np.where(np.abs(arr) > ep)[0]


def sbml_to_reaction_list(sbml_file):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(sbml_file)
    model = document.getModel()

    list_reactions = []

    for i in range(model.getNumReactions()):
        reaction = model.getReaction(i)
        reactants = []
        products = []

        for reactant in reaction.getListOfReactants():
            species_id = reactant.getSpecies()
            stoichiometry = reactant.getStoichiometry()
            is_boundary = model.getSpecies(species_id).getBoundaryCondition()
            if not is_boundary:
                reactants.append([species_id, stoichiometry])

        for modifier in reaction.getListOfModifiers():
            species_id = modifier.getSpecies()
            stoichiometry = 1
            is_boundary = model.getSpecies(species_id).getBoundaryCondition()
            if not is_boundary:
                reactants.append([species_id, stoichiometry])

        if len(reactants) == 0:
            reactants.append(['', 1])

        for product in reaction.getListOfProducts():
            species_id = product.getSpecies()
            stoichiometry = product.getStoichiometry()
            is_boundary = model.getSpecies(species_id).getBoundaryCondition()
            if not is_boundary:
                products.append([species_id, stoichiometry])

        if len(reaction.getListOfProducts()) == 0:
            products.append(['', 1])

        list_reactions.append([reactants, products])

        if reaction.getReversible():
            list_reactions.append([products, reactants])
    
    return list_reactions


def _transpose(matrix): # transpose a list of lists
    return list(map(list, zip(*matrix)))


def _compute_nullspace(mat, mode='dense'):
    if mode == 'dense':
        return null_space(mat)
    else:
        raise ValueError("Invalid mode.")


def _deep_compare(lst1, lst2):
    if type(lst1) != type(lst2):
        print("Type mismatch:", lst1, lst2)
        return False
    if isinstance(lst1, list) and isinstance(lst2, list):
        if len(lst1) != len(lst2):
            print("Length mismatch:", lst1, lst2)
            return False
        for i, (sub_lst1, sub_lst2) in enumerate(zip(lst1, lst2)):
            if not _deep_compare(sub_lst1, sub_lst2):
                print(f"Mismatch found in list at index {i}: {sub_lst1} != {sub_lst2}")
                return False
        return True
    elif isinstance(lst1, np.ndarray) and isinstance(lst2, np.ndarray):
        if lst1.shape != lst2.shape:
            print("Shape mismatch:", lst1.shape, lst2.shape)
            return False
        comparison = np.array_equal(lst1, lst2)
        if not comparison:
            print("Array value mismatch")
        return comparison
    else:
        comparison = lst1 == lst2
        if isinstance(comparison, np.ndarray):
            if not comparison.all():
                print("Direct value mismatch in arrays:", lst1, lst2)
            return comparison.all()
        else:
            if not comparison:
                print("Direct value mismatch:", lst1, lst2)
            return comparison



class ReactionSystem:
    """
    Represents a reaction system. Implements functions for finding RPA properites of a reaction system.

    Attributes:
    - species (list): The list of species.
    - nr (int): The number of reactions.
    - ns (int): The number of species.
    - ind (dict): The dictionary mapping species to their indices.
    - rmat (numpy.ndarray): The reaction matrix.
    - pmat (numpy.ndarray): The product matrix.
    - smat (numpy.ndarray): The stoichiometric matrix.
    - cmat (numpy.ndarray): The nullspace matrix of smat.
    - dmat (numpy.ndarray): The nullspace matrix of smat transposed.

    Methods:
    - enumerate_labeled_buffering_structures: Enumerates labeled buffering structures.
    - enumerate_labeled_buffering_structures_by_name: Enumerates labeled buffering structures (pecies are indicated by names).
    - enumerate_affected_subnetworks: Enumerates subnetworks that are affected by the perturbation of each parameter.
    - compute_influence_index: Computes the influence index of a subnetwork.
    - is_output_complete: Checks if a given subnetwork is output-complete. 
    - find_reactions_to_add: For a given subnetwork, it finds reactions to be added to make the subnetwork output-complete.
    - write_species_to_file: Writes the species to a file.
    - get_emergent_cycles: Returns the emergent cycles of a subnetwork.
    - get_emergent_cycles_symb: Returns the emergent cycles of a subnetwork computed using sympy.
    - get_emergent_conserved_quantities: Returns the emergent conserved quantities of a subnetwork.
    - get_emergent_conserved_quantities_symb: Returns the emergent conserved quantities of a subnetwork computed using sympy.
    - get_lost_conserved_quantities: Returns the lost conserved quantities of a subnetwork.
    - get_lost_conserved_quantities_symb: Returns the lost conserved quantities of a subnetwork computed using sympy.
    - find_integrators_from_lbs: Returns the integrators corresponding to a labeled buffering structure.
    - find_integrators_from_bs: Returns the integrators for a buffering structure. 
    """

    def __init__(self, reactions_input, input_type='string', show_info=False):
        """
        Initializes a ReactionSystem object.

        Args:
            reactions_input (str or list): The input reactions. If input_type is 'string', reactions_input should be 
                a string representing the reactions. If input_type is 'list', reactions_input should be a list of reactions.
                If input_type is 'sbml', reactions_input should be the filename of an SBML file.
            input_type (str): The type of input reactions. Default is 'string'. Other options are 'list' and 'sbml'.
            show_info (bool): Whether to display additional information during initialization. Default is False.

        Returns:
            None

        """

        if isinstance(reactions_input, str) and input_type == 'string':
            reactions = parse_string_to_reaction_list(reactions_input)
            self.nr = len(reactions)
            r_list_c = [ [from_str_to_list(reactions[i][j]) for i in range(self.nr) ] for j in range(2) ] 
        elif isinstance(reactions_input, list): 
            reactions = reactions_input
            self.nr = len(reactions)
            r_list_c = [ [from_str_to_list(reactions[i][j]) for i in range(self.nr) ] for j in range(2) ] 
        elif input_type == 'sbml':
            r_list_c = _transpose(sbml_to_reaction_list(reactions_input))
            self.nr = len(r_list_c[0])
        else:
            raise NotImplementedError("unrecognized file format")

        if show_info:
            print('reactions\n', r_list_c)

        self.species = extract_species(flatten_list(r_list_c))
        self.ns = len(self.species)
        self.ind = { self.species[i] : i for i in range(len(self.species)) }

        if show_info:
            print("\nspecies list")
            print(self.ind)

        self.rmat = np.transpose(np.array([ self._complex_to_vec(r) for r in r_list_c[0]]))
        self.pmat = np.transpose(np.array([ self._complex_to_vec(r) for r in r_list_c[1]]))
        self.smat = self.pmat - self.rmat 
        self.cmat = _compute_nullspace(self.smat).T
        self.dmat = _compute_nullspace(self.smat.T).T
        # row vectors correspond to cycles and conserved quantities

        tolerance = 1e-10
        self.cmat[np.abs(self.cmat) < tolerance] = 0
        self.dmat[np.abs(self.dmat) < tolerance] = 0

        if show_info:
            print("\nnum. of reactioins:", self.nr)
            print("num. of species:", self.ns)
            print("dim ker S =", len(self.cmat))
            print("dim coker S =", len(self.dmat))
        
        self._lbs_list = None
        self._affected_subnetworks = None

        return

    def write_species_to_file(self, filename):
        """
        Write the species to a file.

        Args:
            filename (str): The name of the file to write the species to.

        Returns:
            None
        """
        current_directory = os.getcwd()
        filename = os.path.join(current_directory, filename)
        with open(filename, 'w') as f:
            for s in self.species:
                f.write(s + '\n')
        return

    def _complex_to_vec(self, c):
        # c = [["A",1],["B",2]], [["C",2], ["B",1]]        
        c = [cc for cc in c if cc[0]!=''] 
        nums = [ [self.ind[s[0]], s[1]] for s in c]
        vec = np.zeros(len(self.species), dtype='float64')
        for num in nums: 
            vec[num[0]] = num[1]
        return vec

    def compute_influence_index(self, subn, verbose=False):
        """
        Computes the influence index of a given subnetwork.

        Parameters:
        - subn: A tuple containing the subnetwork's species (xs) and reactions (rs).
        - verbose: A boolean indicating whether to print additional information.

        Returns:
        - int: The influence index of the subnetwork.
        """
        xs, rs = subn

        smat = self.smat        
        xall = list(range(self.ns))
        rall = list(range(self.nr))
        pr = proj_mat(rs, rall)
        pmbar = proj_mat_comp(xs, xall)

        if verbose:
            print('len(rall):', len(rall) )
            print('dimcoker(smat):', dimcoker(smat))
            print('dimker(np.dot(smat, pr)):', dimker(np.dot(smat, pr)))
            print('dimcoker(np.dot(pmbar, smat)):', dimcoker(np.dot(pmbar, smat)))
            print('xs:', xs)
            print('xall:', xall)
            print('rs:', rs)
            print('rall:', rall)
            print('pmbar\n', pmbar)
            print('pmbar @ smat:\n')
            print(pmbar @ smat)

        return len(rall) + dimcoker(smat) - dimker(smat @ pr) - dimcoker(pmbar @ smat)

    def is_output_complete(self, subn):
        """
        Check if a given subnetwork is output-complete.

        Parameters:
        subn (tuple): A tuple containing two lists: xs and rs.
                  xs (list): List of species indices.
                  rs (list): List of reaction indices.

        Returns:
        bool: True if the subnetwork is output-complete, False otherwise.
        """
        xs, rs = subn

        aa = np.array( [ self.rmat[i,:] for i in xs ] )
        r_indices = np.array(flatten_list([ np.nonzero(a)[0] for a in aa]))
        r_indices = np.unique(r_indices)

        return _contains_all(rs, r_indices.tolist())

    def _compute_amat(self, verbose=False):
        drdx = _randomize_nonzero_entries(self.rmat.T)

        cmat = self.cmat
        dmat = self.dmat
        
        ns, nr = self.ns, self.nr
        nc, nd = cmat.shape[0], dmat.shape[0] # indices are: c[alpha, A], d[alpha_bar, i]
        dim = ns + nc
        
        if verbose:
            print(f"nx,nr,nc,nd: {ns, nr, nc, nd}")
            print("drdx:\n", drdx)
            print("cmat:\n", cmat)
            print(cmat.shape)
            print("dmat:\n", dmat)

        assert ns + nc == nr + nd, "Fredholm's theorem violated"
    
        amat = np.zeros((dim, dim))
        for i in range(dim):
            for j in range(dim):
                if i < nr and j < ns:
                    amat[i, j] = drdx[i, j]
                elif i < nr and j >= ns:
                    amat[i, j] = cmat[j - ns, i]
                elif i >= nr and j < ns:
                    amat[i, j] = dmat[i - nr, j]
        
        return amat

    def _compute_response_mat(self, verbose=False):        
        ns = self.ns
        nc = self.cmat.shape[0]
        dim = ns + nc

        amat = self._compute_amat()

        # det = np.linalg.det(amat)
        # print("det A=", det)
        # if np.abs(det) < 1e-10:
        #     raise ValueError("Matrix is not invertible; determinant is close to zero: ", det)

        ainv = np.linalg.inv(amat)

        if verbose:
            print("amat\n", amat)
            print("ainv\n", ainv)
        
        s_resp = ainv[:ns, :]
        r_resp = self.cmat.T @ ainv[ns:dim, :]
        
        return s_resp, r_resp

    def find_reactions_to_add(self, xs, rs):
        """
        Finds the reactions to add to make the given subnetwork output-complete.

        Args:
            xs (list): List of input species.
            rs (list): List of output species.

        Returns:
            numpy.ndarray: Array of indices of reactions to be added. 
        """
        
        if self.is_output_complete([xs, rs]):
            return np.array([], dtype=int)
        
        rmat = self.rmat
        list_ = np.sum(rmat[np.array(xs)], axis=0)
        affected_rs = np.nonzero(list_ > 0)[0]

        complement = list(set(affected_rs) - set(rs))

        return sorted(complement)

    def enumerate_affected_subnetworks(self):
        """
        Enumerates the affected subnetworks by the perturbation of each reaction parameter.

        Parameters:
        - verbose (bool): If True, prints the responded species and reactions.

        Returns:
        - affected_subnetworks (list): A list of affected subnetworks, whose i-th subnetwork is represents subnetworks affected by the perturbation of the i-th reaction parameter.
        """

        if self._affected_subnetworks is not None:
            return self._affected_subnetworks
        ns = self.ns
        nc = self.cmat.shape[0]
        dim = ns + nc

        s_resp, r_resp = self._compute_response_mat()

        responded_species = []
        for j in range(dim):
            non_zero_indices = _nonzero_with_error(s_resp[:, j])
            responded_species.append(non_zero_indices)

        responded_reactions = []
        for j in range(dim):
            non_zero_indices = _nonzero_with_error(r_resp[:, j])
            responded_reactions.append(non_zero_indices)

        affected_subnetworks = []
        for i in range(dim):
            xs = responded_species[i]
            rs = responded_reactions[i]
            affected_subnetworks.append([xs, rs])
        
        self._affected_subnetworks = affected_subnetworks

        return affected_subnetworks

    def enumerate_labeled_buffering_structures(self, ntrial=2, verbose=False):
        """
            Enumerates labeled buffering structures of a given reaction network.

            Parameters:
            - ntrial (int): The number of trials to perform for identifying affected subnetworks. Default is 2.
            - verbose (bool): Whether to print additional information during the enumeration process. Default is False.

            Returns:
            - lbs_list (list): A list of labeled buffering structures, where each structure is represented as a nested list.
        """

        if self._lbs_list is not None:
            return self._lbs_list

        lbs_would_be_list = [ self.enumerate_affected_subnetworks() for i in range(ntrial) ]

        def is_consistent(lbs_would_be_list):
            return all( [ _deep_compare(lbs_would_be_list[0], el) for el in lbs_would_be_list[1:] ])

        if ntrial>1 and not is_consistent(lbs_would_be_list):
            raise ValueError("The identification of affected subnetworks is not consistent across trials.")

        lbs_would_be = _make_hashable(lbs_would_be_list[0])

        if verbose:
            print('would be lbs')
            for l in lbs_would_be:
                print(l)
        
        if verbose:
            print('hashable form of would be lbs')
            print(_make_hashable(lbs_would_be))

        lbs_duplicates_removed = _remove_duplicates(lbs_would_be)

        added_reactions = []
        for i in range(len(lbs_duplicates_removed)):
            xs = lbs_duplicates_removed[i][0]
            rs = lbs_duplicates_removed[i][1]
            added_reactions.append(tuple(self.find_reactions_to_add(xs, rs)))

        lbs_duplicates_removed_2 = []
        for i in range(len(lbs_duplicates_removed)):
            lbs_duplicates_removed_2.append( lbs_duplicates_removed[i] + (added_reactions[i],) )

        lbs = []
        for i in range(len(lbs_duplicates_removed)):
            reactions = tuple(_get_indices_of_element(lbs_would_be, lbs_duplicates_removed[i]))
            lbs.append( (reactions,) + lbs_duplicates_removed_2[i] )

        lbs_list = [[[el for el in inner_tuple] for inner_tuple in outer_tuple] for outer_tuple in lbs]
        self._lbs_list = lbs_list

        return lbs_list

    def enumerate_labeled_buffering_structures_by_name(self):
        """
        Enumerates the labeled buffering structures (species are indicated by names).

        Returns:
            A list labeled buffering structures, in which species are indicated by names.
        """
        lbs_list = self.enumerate_labeled_buffering_structures()
        return [self.lbs_to_name(l) for l in lbs_list]

    def enumerate_buffering_structures(self):
        """
        Enumerates the buffering structures of a given reaction network.

        Returns:
            list: A list of buffering structures.
        """
        lbs_list = self.enumerate_labeled_buffering_structures()
        bs_list = [ from_lbs_to_bs(l) for l in lbs_list ]
        return bs_list

    def lbs_to_name(self, lbs):
        name = {v: k for k, v in self.ind.items()}
        res = [ lbs[0], list(map(lambda x: name[x], lbs[1])), lbs[2], lbs[3] ]
        return res

    def bs_to_name(self, bs):
        name = {v: k for k, v in self.ind.items()}
        res = [ list(map(lambda x: name[x], bs[0])), bs[1] ]
        return res

    def _get_S11(self, subn):
        return self.smat[np.ix_( *subn )]

    def _get_S12(self, subn):
        xs, rs = subn
        smat = self.smat
        rs_comp = list(set(range(self.nr)) - set(rs))
        return smat[np.ix_(xs, rs_comp)]

    def _get_S21(self, subn):
        xs, rs = subn
        smat = self.smat
        xs_comp = list(set(range(self.ns)) - set(xs))
        return smat[np.ix_(xs_comp, rs)]
    
    def _get_S22(self, subn):
        xs, rs = subn
        smat = self.smat
        xs_comp = list(set(range(self.ns)) - set(xs))
        rs_comp = list(set(range(self.nr)) - set(rs))
        return smat[np.ix_(xs_comp, rs_comp)]

    def _get_reorganized_smat(self, subn):
        s11 = self._get_S11(subn)
        s12 = self._get_S12(subn)
        s21 = self._get_S21(subn)
        s22 = self._get_S22(subn)
        return np.vstack( (np.hstack( (s11, s12) ), np.hstack( (s21, s22) )) )

    def _get_S_tilde(self, subn):
        #
        # S_tilde = (S11 S12 S11)
        #           (S21 S22  0 )
        #
        s11 = self._get_S11(subn)
        n_pad = self.ns - s11.shape[0]
        s11_padded = np.pad(s11, ((0, n_pad), (0, 0)))
        _smat = self._get_reorganized_smat(subn)
        return np.hstack( (_smat, s11_padded) )

    def get_emergent_cycles(self, subn):
        """
        Computes the emergent cycles of a given subnetwork.

        Parameters:
        - subn: A subnetwrok represented as a tuple containing a list of species indices and a list of reaction indices.

        Returns:
        - emergent_cycles: An array containing the emergent cycles of the subnetwork.
        """
        s11 = self._get_S11(subn)
        s21 = self._get_S21(subn)

        if s11.shape[0] == 0:
            kerS11 = np.identity(s11.shape[1])
        else:
            kerS11 = null_space(s11).transpose()
        if kerS11.shape[0] == 0:
            return np.array([])

        if s21.shape[0] == 0:
            kerS21 = np.identity(s21.shape[1])
        else:
            kerS21 = null_space(s21).transpose()
        if kerS21.shape[0] == 0:
            return kerS11
        
        kerS21_orth = get_orthogonal_complement(kerS21)
        return get_intersection(kerS11, kerS21_orth)

    def get_emergent_cycles_symb(self, subn):
        """
        Computes the emergent cycles of a given subnetwork using sympy.

        Parameters:
        - subn: A subnetwrok represented as a tuple containing a list of species indices and a list of reaction indices.

        Returns:
        - emergent_cycles: An array containing the emergent cycles of the subnetwork.
        """
        s11 = sp.Matrix(self._get_S11(subn))
        s21 = sp.Matrix(self._get_S21(subn))

        if s11.rows == 0:
            kerS11 = sp.eye(s11.cols)
        else:
            kerS11 = s11.nullspace()
            kerS11 = sp.Matrix.hstack(*kerS11).transpose()
        if not kerS11:
            return sp.zeros(0, len(subn[1]))

        if s21.rows == 0:
            kerS21 = sp.eye(s21.cols)
        else:
            kerS21 = s21.nullspace()
            kerS21 = sp.Matrix.hstack(*kerS21).transpose()
        if not kerS21:
            return kerS11

        kerS21_orth = get_orthogonal_complement_symb(kerS21)

        return get_intersection_symb(kerS11, kerS21_orth)

    def get_emergent_conserved_quantities(self, subn):
        """
        Computes the emergent conserved quantities for a given subnetwork.

        Args:
        - subn: A subnetwrok represented as a tuple containing a list of species indices and a list of reaction indices.

        Returns:
        - emergent conserved quantities: An array containing the emergent conserved quantities of the subnetwork.
        """
        xs, rs = subn
        S11 = self._get_S11(subn)

        if S11.shape[1] == 0:
            coker_S11 = np.identity(len(xs))
        else:
            coker_S11 = null_space(self._get_S11(subn).T).T

        S_tilde = self._get_S_tilde(subn)
        coker_S_tilde = null_space(S_tilde.T).T
        d1 = coker_S_tilde[:, :len(xs)]

        return get_intersection(coker_S11, get_orthogonal_complement(d1))
    
    def get_emergent_conserved_quantities_symb(self, subn):
        """
        Computes the emergent conserved quantities for a given subnetwork using sympy.

        Args:
        - subn: A subnetwrok represented as a tuple containing a list of species indices and a list of reaction indices.

        Returns:
        - emergent conserved quantities: An array containing the emergent conserved quantities of the subnetwork.
        """
        xs, rs = subn

        S11_np = self._get_S11(subn)
        S11 = sp.Matrix(S11_np.astype(int))

        if S11.cols == 0:
            coker_S11 = sp.eye(len(xs))
        else:
            coker_S11 = S11.T.nullspace()
            coker_S11 = sp.Matrix.hstack(*coker_S11).transpose()
            if coker_S11.rows == 0:
                coker_S11 = sp.zeros(0, len(xs))

        S_tilde_np = self._get_S_tilde(subn)
        S_tilde = sp.Matrix(S_tilde_np.astype(int))

        coker_S_tilde = S_tilde.T.nullspace()

        if len(coker_S_tilde) == 0:
            coker_S_tilde = sp.zeros(0, self.ns)
        else:
            coker_S_tilde = sp.Matrix.hstack(*coker_S_tilde).transpose()

        d1 = coker_S_tilde[:, :len(xs)]

        return get_intersection_symb(coker_S11, get_orthogonal_complement_symb(d1))

    def get_lost_conserved_quantities(self, subn):
        """
        Computes the lost conserved quantities for a given subnetwork.

        Args:
        - subn: A subnetwrok represented as a tuple containing a list of species indices and a list of reaction indices.

        Returns:
        - lost conserved quantities: An array containing the lost conserved quantities of the subnetwork.
        """
        # D_l = coker S \cap coker S_tilde

        S_mat = self._get_reorganized_smat(subn)
        coker_S = _compute_nullspace(S_mat.T).T
        
        S_tilde = self._get_S_tilde(subn)
        coker_S_tilde = null_space(S_tilde.T).T

        return get_intersection(coker_S, get_orthogonal_complement(coker_S_tilde))

    def get_lost_conserved_quantities_symb(self, subn):
        """
        Computes the lost conserved quantities for a given subnetwork using sympy.

        Args:
        - subn: A subnetwrok represented as a tuple containing a list of species indices and a list of reaction indices.

        Returns:
        - lost conserved quantities: An array containing the lost conserved quantities of the subnetwork.
        """
#        S_mat = sp.Matrix(self.smat)
        S_mat = sp.Matrix(self._get_reorganized_smat(subn))

        coker_S = S_mat.T.nullspace()
        coker_S = sp.Matrix.hstack(*coker_S).transpose()

        S_tilde_np = self._get_S_tilde(subn)
        S_tilde = sp.Matrix(S_tilde_np.astype(int))

        coker_S_tilde = S_tilde.T.nullspace()
        coker_S_tilde = sp.Matrix.hstack(*coker_S_tilde).transpose()

        if coker_S_tilde.rows == 0:
            coker_S_tilde = sp.zeros(0, coker_S.cols)

        return get_intersection_symb(coker_S, get_orthogonal_complement_symb(coker_S_tilde))

    def _compute_integrator_matrices(self, subn, pinv_mode='sympy'):

        s11 = sp.Matrix(self._get_S11(subn).astype(int))
        s21 = sp.Matrix(self._get_S21(subn).astype(int))
        s12 = sp.Matrix(self._get_S12(subn).astype(int))
        s22 = sp.Matrix(self._get_S22(subn).astype(int))

        if pinv_mode == 'sympy':
            s11_plus = s11.pinv()
            s21_plus = s21.pinv()
        elif pinv_mode == 'numpy':
            if s11.rows == 0 or s11.cols == 0:
                s11_plus = sp.zeros(s11.cols, s11.rows)
            else:
                s11_plus = sp.Matrix( np.linalg.pinv(np.array(s11.evalf().tolist(), dtype=np.float64)) )
            
            if s21.rows == 0 or s21.cols == 0:
                s21_plus = sp.zeros(s21.cols, s21.rows)
            else:
                s21_plus = sp.Matrix( np.linalg.pinv(np.array(s21.evalf().tolist(), dtype=np.float64)) )
        else:
            raise ValueError("Invalid pinv mode")

        s_p = s22 - s21 * s11_plus * s12   # Generalized Schur complement

        emergent_cycles = self.get_emergent_cycles_symb(subn)
        emergent_conserved_quantities = self.get_emergent_conserved_quantities_symb(subn)

        s21c11 = s21 * emergent_cycles.transpose()

        p_mat = sp.Matrix.hstack( *(s21c11.transpose().nullspace()) ).transpose()
        q_mat = sp.Matrix.hstack(- s21 * s11_plus, sp.eye(s22.rows))

        kerS11 = sp.Matrix.hstack(*s11.nullspace()).transpose() 
        if kerS11.shape[0] == 0:
            kerS11 = sp.zeros(0, s11.cols)

        kerS21 = sp.Matrix.hstack(*s21.nullspace()).transpose() 
        if kerS21.shape[0] == 0:
            kerS21 = sp.zeros(0, s21.cols)

        return [ p_mat * q_mat, 
                 p_mat * s_p, 
                 emergent_conserved_quantities, 
                 emergent_conserved_quantities * s12, 
                 s11_plus, 
                 s11_plus * s12, 
                 kerS11, 
                 s21_plus, 
                 kerS21, 
                 s22
                ]

    def _create_species_str(self, symbol_mode='subscript'):
        if symbol_mode == 'subscript':
            symbol_x_str = " ".join( map(lambda s: "x_{" + s + "}", self.species) )
        elif symbol_mode == 'name':
            symbol_x_str = " ".join( self.species )
        elif symbol_mode == 'index':
            symbol_x_str = " ".join( map(lambda s: "x_{" + str(s) + "}", range(self.ns)) )
        else:
            raise ValueError("Invalid symbol mode")
       
        return symbol_x_str

    def find_integrators_from_bs(self, subn, output_style='latex', symbol_mode='subscript', pinv_mode='sympy'):
        """
        Returns the integrators for a given buffering structure.

        Parameters:
        - subn: A buffering structure, which is a subnetwork with zero influence index.
        - output_style: The style of the output. Valid options are 'latex' and 'symbolic'.
        - symbol_mode: The mode for representing symbols. Valid options are 'subscript', 'name', and 'index'.
        - pinv_mode: The mode for computing the pseudo-inverse. Valid options are 'sympy' and 'numpy'.

        Returns:
        - list: A list of four elements representing the integrators. The type of the elements depends on the parameter `output_style`.
        If we denote the returned list as X, then the integrators are given by:

                $\frac{d}{dt} X[0] = X[1]$

                $\frac{d}{dt} X[2] = X[3]$
        
        Warnings:
        - Using 'sympy' for the pseudo-inverse computation can be slow for large subnetworks. 
        For larger matrices or subnetworks, consider using 'numpy' to significantly reduce computation time.
        - When 'numpy' is selected for `pinv_mode`, the computation uses numerical methods, which 
        might introduce floating-point precision inaccuracies.
        """
        p_q, p_sp, d, c2, s11plus, s11p_s12, kerS11, s21plus, kerS21, s22= self._compute_integrator_matrices(subn, pinv_mode)
        x1, x2, r1, r2, r1_indices = self._prepare_variables(symbol_mode, subn)

        expressions1 = [ d * x1, c2 * r2, p_q * (x1.col_join(x2)), p_sp * r2 ] 

        if output_style == 'latex':
            data = [ (lambda x:sp.latex(sp.simplify(x)))(x) for x in expressions1]
        elif output_style == 'symbolic':
            data = [ sp.simplify(x) for x in expressions1]
        else:
            raise ValueError("Invalid output style")

        return data

    def find_integrators_from_lbs(self, lbs, output_style='latex', symbol_mode='subscript', pinv_mode='sympy'):
        """
        Returns the integrators for a given labeled buffering structure.

        Parameters:
        - lbs (list): A labeled buffering structure.
        - output_style (str, optional): The output style for the integrators. Default is 'latex' If 'symbol' is chosen, it returns sympy expressions.
        - symbol_mode (str, optional): The mode for representing symbols. Default is 'subscript'.
        - pinv_mode: The mode for computing the pseudo-inverse. Valid options are 'sympy' and 'numpy'.

        Returns:
        - list: A list of eights elements representing the integrators. The type of the elements depends on the parameter `output_style`.
        If we denote the returned list as X, then the integrators are given by:

                $\frac{d}{dt} X[0] = X[1]$

                $\frac{d}{dt} X[2] = X[3]$

                $\frac{d}{dt} X[4] = X[5]$

                $\frac{d}{dt} X[6] = X[7]$
        
        Warnings:
        - Using 'sympy' for the pseudo-inverse computation can be slow for large subnetworks. 
        For larger matrices or subnetworks, consider using 'numpy' to significantly reduce computation time.
        - When 'numpy' is selected for `pinv_mode`, the computation uses numerical methods, which 
        might introduce floating-point precision inaccuracies.
        """
        subn = from_lbs_to_bs(lbs)
        p_q, p_sp, d, c2, s11plus, s11plus_s12, kerS11, s21plus, kerS21, s22 = self._compute_integrator_matrices(subn, pinv_mode)
        x1, x2, r1, r2, r1_indices = self._prepare_variables(symbol_mode, subn)
        
        expressions1 = [ d * x1, c2 * r2, p_q * (x1.col_join(x2)), p_sp * r2 ] 

        if output_style == 'latex':
            data1 = [ (lambda x:sp.latex(sp.simplify(x)))(x) for x in expressions1]
        elif output_style == 'symbolic':
            data1 = expressions1
        else:
            raise ValueError("Invalid output style")

        if len(lbs[3]) == 0:  # when \mathcal E_\gamma is empty 
            data2 = [ sp.zeros(0, 1) for _ in range(4)]
            if output_style == 'latex':
                data2 = [ sp.latex(x) for x in data2 ]
                return data1 + data2
            elif output_style == 'symbolic':
                return data1 + data2
            else:
                raise ValueError("Invalid output style")

        # Make a projection matrix to \mathcal E_\gamma
        lbs3_indices = [ r1_indices.index(val) for val in lbs[3] ]
        proj_hat = sp.Matrix.hstack(*[sp.Matrix([0]*i + [1] + [0]*(len(subn[1])-i-1)) for i in lbs3_indices]).transpose()

        kerS11_perp = get_orthogonal_complement_symb(kerS11)
        kerS21_perp = get_orthogonal_complement_symb(kerS21)

        proj_hat_1 = get_intersection_symb(proj_hat, kerS11_perp)
        proj_hat_2 = get_intersection_symb(proj_hat, kerS21_perp)

        expressions2 = [ proj_hat_1 * (s11plus * x1), 
                        proj_hat_1 * (r1 + s11plus_s12 * r2), 
                        proj_hat_2 * (s21plus * x2), 
                        proj_hat_2 * (r1 + s21plus * s22 * r2)
                        ]

        data2 = [ sp.simplify(x) for x in expressions2]

        if output_style == 'latex':
            data2= [ sp.latex(x) for x in data2]
            return data1 + data2
        elif output_style == 'symbolic':
            return data1 + data2
        else:
            raise ValueError("Invalid output style")

    def _prepare_variables(self, symbol_mode, subn):
        symbol_x_str = self._create_species_str(symbol_mode)
        x = sp.Matrix( sp.symbols(symbol_x_str) )

        symbol_r_str = " ".join( map(lambda s: "r_{" + str(s+1) + "}", range(self.nr)) )
        r = sp.Matrix( sp.symbols(symbol_r_str) )

        x1 = sp.Matrix([x[i] for i in subn[0]])
        xs_comp = list(set(range(self.ns)) - set(subn[0]))
        x2 = sp.Matrix([x[i] for i in xs_comp])
        if x2.shape[0] == 0:
            x2 = sp.zeros(0, 1)

        r1 = sp.Matrix([r[i] for i in subn[1]])
        r1_indices = [i for i in range(self.nr) if i in subn[1]]

        rs_comp = list(set(range(self.nr)) - set(subn[1]))
        r2 = sp.Matrix([r[i] for i in rs_comp])
        if r2.shape[0] == 0:
            r2 = sp.zeros(0, 1)

        return x1, x2, r1, r2, r1_indices

    def _create_matrix_affected(self):
        # create a num_reactions by num_species boolean matrix
        matrix_affected_species = np.zeros((self.ns, self.nr), dtype=bool)
        matrix_affected_reactions = np.zeros((self.nr, self.nr), dtype=bool)
        affected_subnetworks = self.enumerate_affected_subnetworks()
        for i in range(self.nr):
            xs, rs = affected_subnetworks[i]
            for x in xs:
                matrix_affected_species[x, i] = True
            for r in rs:
                matrix_affected_reactions[r, i] = True

        return matrix_affected_species, matrix_affected_reactions

    def draw_affected_matrix(self, output_style='show', filename_prefix='plot'):
        """
        Draws and displays or saves the affected species and affected reactions matrices.

        Parameters:
            output_style (str): The output style. Options are 'show' to display the plot or 'pdf' to save the plot as a PDF file. Default is 'show'.
            filename_prefix (str): The prefix for the filename when saving the plot as a PDF file. Default is 'plot'.

        Raises:
            ValueError: If the output_style is not 'show' or 'pdf'.

        Returns:
            None
        """
        matrix_affected_species, matrix_affected_reactions = self._create_matrix_affected()
        sorted_matrix, sorted_row_indices, sorted_col_indices = _sort_rows_and_columns_by_true_count(matrix_affected_species)
        matrix_affected = sorted_matrix
        # print matrix size
        factor = 0.5
        plt.figure(figsize=(factor * self.ns, factor * self.nr))
        plt.title('Affected Species Matrix')
        plt.gca().set_xticks(np.arange(self.nr) + factor, ['' for _ in np.arange(self.nr)], minor=False)
        plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
        # Extract the locations of minor ticks as a list
        minor_tick_locations = plt.gca().xaxis.get_minorticklocs()
        #  remove first and last component of minor_tick_locations
        minor_tick_locations = minor_tick_locations[1:-1]
        # set the minor ticks
        plt.gca().set_xticks(minor_tick_locations, [f'$R_{{{i + 1}}}$' for i in sorted_col_indices], minor=True)
        plt.imshow(matrix_affected, cmap='binary', interpolation='none')
        plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
        # put x ticks on the top
        plt.gca().xaxis.tick_top()

        formatted_species_list = [ f"${s}$" for s in self.species[sorted_row_indices] ]

        plt.gca().set_yticks(np.arange(self.ns) + factor, ['' for _ in np.arange(self.ns)], minor=False)
        plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
        # Extract the locations of minor ticks as a list
        minor_tick_locations = plt.gca().yaxis.get_minorticklocs()
        #  remove first and last component of minor_tick_locations
        minor_tick_locations = minor_tick_locations[2:-1]
        plt.gca().set_yticks(minor_tick_locations, formatted_species_list, minor=True)

        if output_style == 'show':
            plt.show()
        elif output_style == 'pdf':
            filename = f"{filename_prefix}_species.pdf"
            plt.savefig(filename, dpi=300, bbox_inches='tight')
        else:
            raise ValueError("Invalid output_style. Choose 'show' or 'pdf'.")

        sorted_matrix = matrix_affected_reactions[sorted_col_indices][:, sorted_col_indices]
        # start a new figure
        plt.figure(figsize=(factor * self.nr, factor * self.nr))
        # add title to the figure
        plt.title('Affected Reactions Matrix')
        plt.gca().set_xticks(np.arange(self.nr) + factor, ['' for _ in np.arange(self.nr)], minor=False)
        plt.gca().xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
        # Extract the locations of minor ticks as a list
        minor_tick_locations = plt.gca().xaxis.get_minorticklocs()
        #  remove first and last component of minor_tick_locations
        minor_tick_locations = minor_tick_locations[1:-1]
        # set the minor ticks
        plt.gca().set_xticks(minor_tick_locations, [f'$R_{{{i + 1}}}$' for i in sorted_col_indices], minor=True)
        plt.imshow(matrix_affected_reactions, cmap='binary', interpolation='none')
        plt.grid(which='major', color='gray', linestyle='-', linewidth=0.5)
        # put x ticks on the top
        plt.gca().xaxis.tick_top()
        plt.gca().set_yticks(np.arange(self.nr) + factor, ['' for _ in np.arange(self.nr)], minor=False)
        plt.gca().yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
        # Extract the locations of minor ticks as a list
        minor_tick_locations = plt.gca().yaxis.get_minorticklocs()
        #  remove first and last component of minor_tick_locations
        minor_tick_locations = minor_tick_locations[2:-1]
        plt.gca().set_yticks(minor_tick_locations, [f'$R_{{{i + 1}}}$' for i in sorted_col_indices], minor=True)

        if output_style == 'show':
            plt.show()
        elif output_style == 'pdf':
            filename = f"{filename_prefix}_reactions.pdf"
            plt.savefig(filename, dpi=300, bbox_inches='tight')
        else:
            raise ValueError("Invalid output_style. Choose 'show' or 'pdf'.")

        return
