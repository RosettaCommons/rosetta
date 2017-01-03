Tests FilterReportAsPoseExtraScoresMover

Steven Lewis, Cyrus Biotechnology, smlewi@gmail.com

The Mover's purpose is to add constraints to a pose based on its current conformation - in other words constraints that stabilize whatever it looks like now, to reduce movement after changes.
Here we are just looking at what gets printed by the mover with various options in the XML.
This test is more concerned with checking what constraints get made under various option regimes; the add_constraints_to_current_conformation test instead checks how the constraints perform in battle.
