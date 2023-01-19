class Protein:
    '''
    This is a protein object for creating protein-protein interaction network
    '''

    protein_interaction_weighted_score = 0

    def __init__(self, accession: str, interactions: str, score: float = 0) -> 'Protein':
        '''
        Initialization from protein-protein interaction data
        :param accession:
        :param interactions:
        '''
        self.accession = accession
        if interactions != None:
            self.interactions = interactions.split('|')
        else:
            self.interactions = None
        self.score = score



    def __str__(self):
        return ': '.join(['Accession', self.accession, 'Score', str(self.score), 'Interactions', '|'.join(self.interactions)])



class ProteinInteraction:
    '''
    The ProteinInteraction class control how information between proteins in the interaction network flows
    '''

    proteins = {}

