# demo_clean_dummy_data.py
"""
一个极简 demo：
用几条假数据跑一遍 clean_reactions，
打印清洗前后反应条目数和部分示例。
"""

from chemrxn_cleaner import clean_reactions, has_product_filter, ReactionRecord


def main():
    # 构造几条简单的 reaction SMILES（有对的有错的）
    raw_reactions = [
        # valid-ish
        "CCO.CCBr>NaOH>CCOC",
        "CC(=O)Cl.NH3>>CC(=O)NH2",
        "c1ccccc1Br>Mg>c1ccccc1MgBr",
        # invalid format
        "CCO>>",  # 没有 product
        ">>>",  # 空反应
        "this_is_not_smiles>NaOH>CCOC",
    ]

    print(f"Raw reactions: {len(raw_reactions)}")

    cleaned: list[ReactionRecord] = clean_reactions(
        raw_reactions,
        filters=[has_product_filter],  # 这里只用一个最简单的 filter 做示意
    )

    print(f"Cleaned reactions: {len(cleaned)}\n")

    print("Examples after cleaning:")
    for rec in cleaned:
        print("-" * 40)
        print(f"Raw:       {rec.raw}")
        print(f"Reactants: {rec.reactants}")
        print(f"Reagents:  {rec.reagents}")
        print(f"Products:  {rec.products}")


if __name__ == "__main__":
    main()
