from __future__ import print_function
import trueskill

# The output of this program should match the output of the TrueSkill
# calculator at:
#
#   http://atom.research.microsoft.com/trueskill/rankcalculator.aspx
#
# (Select game mode "custom", create 4 players each on their own team,
# check the second "Draw?" box to indicate a tie for second place,
# then click "Recalculate Skill Level Distribution".  The mu and sigma
# values in the "after game" section should match what this program
# prints.

# The objects we pass to AdjustPlayers can be anything with skill and
# rank attributes.  We'll create a simple Player class that has
# nothing else.

class Player(object):
  pass

# Create four players.  Assign each of them the default skill.  The
# player ranking (their "level") is mu-3*sigma, so the default skill
# value corresponds to a level of 0.

alice = Player()
alice.skill = (25.0, 25.0/3.0)

bob = Player()
bob.skill = (25.0, 25.0/3.0)

chris = Player()
chris.skill = (25.0, 25.0/3.0)

darren = Player()
darren.skill = (25.0, 25.0/3.0)

# The four players play a game.  Alice wins, Bob and Chris tie for
# second, Darren comes in last.  The actual numerical values of the
# ranks don't matter, they could be (1, 2, 2, 4) or (1, 2, 2, 3) or
# (23, 45, 45, 67).  All that matters is that a smaller rank beats a
# larger one, and equal ranks indicate draws.

alice.rank = 1
bob.rank = 2
chris.rank = 2
darren.rank = 4

# Do the computation to find each player's new skill estimate.

trueskill.AdjustPlayers([alice, bob, chris, darren])

# Print the results.

print(" Alice: mu={0[0]:.3f}  sigma={0[1]:.3f}".format(alice.skill))
print("   Bob: mu={0[0]:.3f}  sigma={0[1]:.3f}".format(bob.skill))
print(" Chris: mu={0[0]:.3f}  sigma={0[1]:.3f}".format(chris.skill))
print("Darren: mu={0[0]:.3f}  sigma={0[1]:.3f}".format(darren.skill))
