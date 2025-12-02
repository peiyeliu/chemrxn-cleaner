# chemrxn_cleaner/ml/dataset.py
from typing import List, Dict, Any
from chemrxn_cleaner.types import ReactionRecord
from torch.utils.data import Dataset


class ForwardReactionDataset(Dataset):
    """
    Minimal forward-prediction dataset built on top of ReactionRecord.
    """

    def __init__(
        self,
        records: List[ReactionRecord],
        use_agents: bool = True,
    ):
        self.records = records
        self.use_agents = use_agents

    def __len__(self) -> int:
        return len(self.records)

    def _make_input_smiles(self, r: ReactionRecord) -> str:
        # left side = reactants [+ (reagents/solvents/agents)]
        left_parts = list(r.reactant_smiles)
        if self.use_agents:
            left_parts += r.reagent_smiles
        return ".".join(left_parts)

    def __getitem__(self, idx: int) -> Dict[str, Any]:
        r = self.records[idx]
        x = self._make_input_smiles(r)
        y = ".".join(r.product_smiles)
        return {
            "input_smiles": x,
            "target_smiles": y,
            "reaction_id": r.reaction_id,
            "meta": {
                "temperature_c": r.temperature_c,
                "time_hours": r.time_hours,
                "solvents": r.solvents,
                "catalysts": r.catalysts,
                "bases": r.bases,
                "additives": r.additives,
                "source": r.source,
            },
        }
