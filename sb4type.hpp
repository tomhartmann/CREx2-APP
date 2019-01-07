#ifndef _SB4TYPE_HPP_
#define _SB4TYPE_HPP_

#include <string>
#include <vector>
#include <iterator>
#include <set>

#include "genom.hpp"
#include "rearrangements.hpp"
using namespace std;

/*
 * multiplication chromosomes of genomes
 * @param[in]: pi: first factor
 * @param[in]: sigma: second factor
 * @return: genome that contains as chromosome the product of pi°sigma; nmap and circ are inherited from sigma
 */
genom multiply(genom pi, genom sigma);

/*
 * calculates a sorting scenario to obtain target from start using only 4-type rrrmts
 * @param[in] start: starting genome of sorting scenario
 * @param[in] target: target genome of sorting scenario
 * @param[in/out] sceanrio: ordered sceario of 4-type rearrangment operations
 */
void sort_by_4type(genom start, genom target, rrrmt* &scenario);

/*
 * calculates a sorting scenario to obtain a permutation g from the identity permutation by only using 4type rrrmts
 * @param[in] g: target genome
 * @param[in/out] rrrmts: rrrmts that sort the genomes of scen
 */
void sort_by_4type(genom g, vector<rrrmt* > &rrrmts);

/*
 * applies all inversions, transpositions, inverse transpositions, and TDRLs to genome g
 * 		and checks whether or not the sorting algorithm uniquely finds the rrrmt
 * @param[in] genom g: just a genome that can be rearranged
 *
 * 	g must still be iota!!!
 */
void test_sorting_algorithm(genom g);

/*
 * applies all inversions, transpositions, inverse transpositions, and TDRLs to genome g
 * 		and checks whether or not the sorting algorithm uniquely finds the rrrmt
 * @param[in] genom g: just a genome that can be rearranged
 *
 *	here g can be an arbirary genome
 */
void test_sorting_algorithm_pairwise(genom g);

/*
 * merges two sets first and second increasingly such that every element of second is inverted
 * @param[in/out] result: the resulting set containing all elements
 * @param[in] first: first set of integer
 * @para[in] second: second set of integer
 */
void oplus(set<int> &result, set<int> first, set<int> second);




#endif// _SB4TYPE_HPP_
