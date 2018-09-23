# EntropyFilter

The script applies several filters to the data. First, it filters out those species which have a percentage of undefined (?) characters above a certain cut-off. Next, it selects those species which are over-lumped. The names of these species should be listed in a separate txt file. The third filtering step is the most important, and is essential to the whole method. This filtering step involved calculating entropy for each character. A column of character values is extracted from the double-filtered data set. Those characters are filtered out, which contain a certain percentage of undefined characters, just as with the row filtering criterion. Shannon entropy is calculated for each of the characters, minus the undefined states of a given character. Mixed characters, such as {0,1} are treated as separate characters (thus, 0, 1, and {0,1} count as three states of a given character). Shannon entropy is calculated in the following manner for a given character:

H=∑_(i=0)^n▒〖〖-p〗_i 〖log〗_n p_i 〗

Where n is equal to the number of character states, pi is equal to the probability of observing state i of the given character, and is equal to (the number of occurrences of state i)/(the total number of occurrences for a given character).
