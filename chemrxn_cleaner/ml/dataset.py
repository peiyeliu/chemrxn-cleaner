# chemrxn_cleaner/ml/dataset.py
from typing import Any, Dict, List

from torch.utils.data import Dataset

from chemrxn_cleaner.types import ReactionRecord


class ForwardReactionDataset(Dataset):
    """Minimal forward-prediction dataset built on top of ReactionRecord."""

    def __init__(
        self,
        records: List[ReactionRecord],
        use_agents: bool = True,
    ):
        """Initialize the dataset.

        Args:
            records: Reaction records providing inputs and targets.
            use_agents: When True, include reagents/agents on the input side.
        """
        self.records = records
        self.use_agents = use_agents

    def __len__(self) -> int:
        """Return the number of reactions in the dataset."""
        return len(self.records)

    def _make_input_smiles(self, r: ReactionRecord) -> str:
        """Assemble the left-hand-side SMILES string for a record.

        Args:
            r: Reaction record providing input SMILES pieces.

        Returns:
            Dot-joined SMILES string for model input.
        """
        # left side = reactants [+ (reagents/solvents/agents)]
        left_parts = list(r.reactant_smiles)
        if self.use_agents:
            left_parts += r.reagent_smiles
        return ".".join(left_parts)

    def __getitem__(self, idx: int) -> Dict[str, Any]:
        """Return model-ready input/target pair and metadata for an index.

        Args:
            idx: Integer index into the records list.

        Returns:
            Dictionary containing ``input_smiles``, ``target_smiles``, and
            metadata useful for downstream models.
        """
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
