package dbOptimizer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class Optimizer {
	final int r;// = 1;
	final int t;// = 2;
	final int l;// = 1;
	final int m;// = 16;
	final int a;// = 2;
	final int f;// = 4;
	final StringBuilder output = new StringBuilder(
			"==================================================================\n");
	final List<double[]> queries = new LinkedList<double[]>();

	public Optimizer(final String queryFile, final String configFile) {
		final Map<String, Integer> configs = new HashMap<String, Integer>();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(configFile));
			String line = reader.readLine();
			while (line != null && !line.isEmpty()) {
				final String[] splits = line.split("=");
				configs.put(splits[0].trim(), Integer.valueOf(splits[1].trim()));
				line = reader.readLine();
			}
			r = configs.get("r");
			t = configs.get("t");
			l = configs.get("l");
			m = configs.get("m");
			a = configs.get("a");
			f = configs.get("f");
			reader.close();
			reader = new BufferedReader(new FileReader(queryFile));
			line = reader.readLine();
			while (line != null) {
				final String[] splits = line.split(" +");
				final double[] oneQuery = new double[splits.length];
				for (int i = 0; i < splits.length; i++) {
					oneQuery[i] = Double.valueOf(splits[i]);
				}
				queries.add(oneQuery);
				line = reader.readLine();
			}
		} catch (final Exception e) {
			throw new RuntimeException(e);
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		// System.out.println("r= " + r + ", t= " + t + ", l= " + l + ", m= " +
		// m
		// + ", a= " + a + ", f= " + f);
	}

	void optimize() {
		for (final double[] selectivity : queries) {
			Combo[] arr = generateSubsets(selectivity, selectivity.length);
			searchOptimal(arr);
			output(arr[arr.length - 1], selectivity);
		}
		System.out.println(output.toString());
	}

	private Combo[] generateSubsets(double[] selectivity, int length) {
		ArrayList<Combo> allSubset = new ArrayList<Combo>();
		// find the stopping point for of iterating all subsets
		int stopMask = 1;
		int stop = 0; // the final subset represent by "11111..1"
		for (int i = 0; i < length; i++) {
			stop += stopMask;
			stopMask <<= 1;
		}
		// generate all initial subset
		int iter = 1;
		while (iter != stop + 1) {
			Combo curr = new Combo();
			curr.index = iter;
			// total selectivity and number of terms
			int numSelect = 0;
			double totalSelect = 1;
			int mask = 1;
			for (int i = 0; i < length; i++) {
				if ((iter & mask) != 0) {
					numSelect++;
					totalSelect *= selectivity[i];
				}
				mask <<= 1;
			}
			curr.select = totalSelect;
			curr.n = numSelect;
			// store c-metric to the leftmostC, and add d-metric to the list
			double p = curr.select;
			double fCost = computeFcost(curr);
			curr.leftmostC = new Pair((p - 1) / fCost, p);
			curr.leastD = new ArrayList<Pair>();
			curr.leastD.add(new Pair(fCost, p));
			// the cost of all & plan
			double cost = numSelect * (r + f) + (numSelect - 1) * l + t + m
					* (curr.select <= 0.5 ? curr.select : 1 - curr.select)
					+ curr.select * a;
			// System.out.println("& Cost " + cost);
			double nonBranchCost = numSelect * (r + f) + (numSelect - 1) * l
					+ a;
			// System.out.println("NonB Cost " + nonBranchCost);
			if (nonBranchCost < cost) {
				cost = nonBranchCost;
				curr.nonBranch = true;
			}
			curr.cost = cost;
			allSubset.add(curr);
			// System.out.println(curr);
			iter++; // to the next subset
		}
		Combo[] result = new Combo[allSubset.size()];
		return allSubset.toArray(result);
	}

	private void searchOptimal(Combo[] subsets) {
		// System.out.println(subsets.length);
		Combo[] optimal = Arrays.copyOf(subsets, subsets.length); // store the
																	// optimal
																	// solutions
																	// from dp
																	// algorithm
		for (int i = 1; i < optimal.length + 1; i++) {
			for (int j = 1; j < subsets.length + 1; j++) {
				if (j != i && (i & j) == 0) {
					// System.out.println("left: " + j + " right: " + i);
					Combo left = subsets[j - 1];
					Combo right = optimal[i - 1];
					// check if the plan is suboptimal
					if (isSuboptimal(left, right)) {
						continue;
					}
					double leftF = computeFcost(left);
					// this plan could be an optimal plan
					// compute the new cost
					double newCost = 0; // newCost = fcost(left) + mq + pC;
										// p is combined selec of left
					newCost = leftF
							+ m
							* (left.select <= 0.5 ? left.select
									: 1 - left.select) + left.select
							* right.cost;
					Combo union = optimal[i + j - 1];
					// System.out.println(union);
					// System.out.println("New cost is " + newCost);
					if (newCost < union.cost) {
						// System.out.println("Find a better solution");
						// update the optimal plan
						union.cost = newCost;
						union.left = left;
						union.right = right;
						union.leftmostC = left.leftmostC; // the leftmost
															// &-term is the
															// left term
						// update the list of best d-metrics
						Pair leftD = new Pair(leftF, left.select);
						union.leastD.add(leftD);

						// LinkedList<Pair> newD = new LinkedList<Pair>();
						// boolean addLeft = true;
						// for (Pair d : right.leastD) {
						// if (d.compareTo(leftD) < 0) {
						// addLeft = false;
						// } else if (d.compareTo(leftD) > 0) {
						// continue;
						// }
						// newD.add(d);
						// }
						// if (addLeft) {
						// newD.add(leftD);
						// }
					}
				}
			}
		}
	}

	/**
	 * Check if the plan is suboptimal by two ways
	 * 
	 * @param left
	 *            : the left &-term that attached to left side of &&
	 * @param right
	 *            : the right side is the current optimal plan
	 * @return
	 */

	private boolean isSuboptimal(Combo left, Combo right) {
		double fLeft = computeFcost(left);
		double p = left.select;
		// compare c-metrics
		Pair cLeft = new Pair((p - 1) / fLeft, p); // c-metric of left &-term
		if (cLeft.compareTo(right.leftmostC) > 0) {
			return true;
		}
		// compare d-metrics
		if (left.select <= 0.5) {
			// check d-metric of L against all possible candidate in R
			Pair dLeft = new Pair(fLeft, p); // d-metric of left &-terms
			for (Pair d : right.leastD) {
				if (dLeft.compareTo(d) > 0) {
					return true;
				}
			}
		}
		return false;
	}

	private double computeFcost(Combo subset) {
		return subset.n * (r + f) + (subset.n - 1) * l + t;
	}

	/**
	 * Print the final optimal plan as a program
	 * 
	 * @param optimal
	 *            : the optimal plan generated by dp
	 * @param select
	 */

	private void output(Combo plan, double[] select) {
		int length = select.length;
		for (final double oneSelect : select) {
			output.append(oneSelect + " ");
		}
		output.append("\n");
		output.append("------------------------------------------------------------------\n");
		// print the branch
		// System.out.println("optimal: " + curr);
		generateC(plan, length);
		output.append("------------------------------------------------------------------\n");
		// print the cost
		output.append("cost: " + String.format("%,.1f", plan.cost) + "\n");
		output.append("==================================================================\n");
	}

	private void generateC(final Combo plan, final int length) {
		output.append("if(");
		final Combo nonBranch = printCondition(plan, length, true);
		output.append("){");
		output.append('\n');
		if (nonBranch != null) {
			output.append("\tanswer[j] = i;\n");
			output.append("\tj+=");
			printCondition(nonBranch, length, false);
			output.append(";\n");
			output.append("}\n");
		} else {
			output.append("\tanswer[j++] = i;\n}\n");
		}
	}

	private Combo printCondition(final Combo plan, final int length,
			final boolean findingNonBranch) {
		if (plan.left == null) {
			if (findingNonBranch && plan.nonBranch) {
				return plan;
			}
			if (plan.n > 1) {
				output.append("(");
			}
			int i = 0;
			int mask = 1;
			boolean firstTerm = true;
			for (; i < length; i++) {
				if ((plan.index & mask) != 0) {
					if (!firstTerm) {
						output.append("&");
					}
					firstTerm = false;
					output.append(printFunction(i + 1));
				}
				mask <<= 1;
			}
			if (plan.n > 1) {
				output.append(")");
			}
			return null;
		} else {
			printCondition(plan.left, length, false);
			if (!findingNonBranch) {
				output.append("&&");
				printCondition(plan.right, length, false);
				return null;
			} else {
				if (plan.right.left == null) {
					if (plan.right.nonBranch) {
						return plan.right;
					}
				}
				output.append("&&");
				return printCondition(plan.right, length, true);
			}
		}
	}

	private String printFunction(int num) {
		StringBuilder sb = new StringBuilder();
		sb.append('t');
		sb.append(num);
		sb.append('[');
		sb.append('o');
		sb.append(num);
		sb.append("[i]");
		sb.append(']');
		return sb.toString();
	}

	/**
	 * c-metric or d-metric provide compareTo method to compare two metrics: 0
	 * stands for not comparable
	 */
	private static class Pair {
		Pair(double score, double p) {
			super();
			this.score = score;
			this.p = p;
		}

		int compareTo(Pair other) {
			if (this.score < other.score && this.p < other.p) {
				return -1;
			} else if (this.score > other.score && this.p > other.p) {
				return 1;

			} else {
				// not comparable
				return 0;
			}
		}

		private double score;
		private double p;
	}

	private static class Combo {
		int index;
		int n;
		// selectivity of all terms
		double select;
		// if non-branch is used
		boolean nonBranch;
		// current cost for the subset
		double cost;
		// the c-metric of leftmost term
		Pair leftmostC;
		List<Pair> leastD;
		Combo left;
		Combo right;

		public Combo() {
			this.index = 0;
			this.n = 0;
			this.select = 0;
			this.nonBranch = false;
			this.cost = 0;
			this.leftmostC = null;
			this.leastD = null;
			this.left = null;
			this.right = null;
		}

		@Override
		public String toString() {
			return "Combo [index=" + index + ", n=" + n + ", selec=" + select
					+ ", nonBranch=" + nonBranch + ", cost=" + cost + "]";
		}
	}

	public static void main(String[] args) {
		new Optimizer(args[0], args[1]).optimize();
	}
}
