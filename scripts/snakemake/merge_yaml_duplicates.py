import yaml
from typing import Any, Dict


def deep_merge_dicts(base: Dict[str, Any], incoming: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge dictionaries. Values from `incoming` merge into `base`.
    - If both values are dicts, merge recursively
    - Otherwise, the `incoming` value overwrites
    """
    for key, value in incoming.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, dict):
            base[key] = deep_merge_dicts(dict(base[key]), value)
        else:
            base[key] = value
    return base


class DuplicateMergeLoader(yaml.SafeLoader):
    """YAML loader that deep-merges duplicate mapping keys instead of last-one-wins."""


def _construct_mapping_with_merge(loader: yaml.Loader, node: yaml.nodes.MappingNode, deep: bool = False):
    # Ensure any YAML merges (<<:) are handled
    loader.flatten_mapping(node)
    mapping: Dict[str, Any] = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        value = loader.construct_object(value_node, deep=deep)
        if key in mapping and isinstance(mapping[key], dict) and isinstance(value, dict):
            mapping[key] = deep_merge_dicts(mapping[key], value)
        else:
            mapping[key] = value
    return mapping


DuplicateMergeLoader.add_constructor(  # type: ignore[attr-defined]
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,  # type: ignore[attr-defined]
    _construct_mapping_with_merge,
)


def load_yaml_with_merged_duplicates(path: str) -> Dict[str, Any]:
    """Load a YAML file and deep-merge duplicate keys.

    Args:
        path: Path to YAML file

    Returns:
        A dict with duplicates deep-merged.
    """
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.load(fh, Loader=DuplicateMergeLoader)
        if data is None:
            return {}
        if not isinstance(data, dict):
            raise TypeError("Top-level YAML must be a mapping (dict)")
        return data 