/*
Milush Yanev
10.2.2019
CS 4200
N-Queen Problem
*/


package nqueenproblem;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Scanner;
import java.util.concurrent.TimeUnit;

/**
 *
 * @author milus
 */


class GeneticAlgorithm {
    public int[] solutionGenerator(int nQeens, int populationSize, double mutation, int numOfGenerations) {
        
        // each one should get a mate.
        populationSize = populationSize - (populationSize % 2); 

        int[][] population = generatePopulation(nQeens, populationSize);

        //set Fitness to max Fitness
        int Fitness = getMaxFitness(nQeens);

        //loop until max num of generations is reached
        for (int i = 0; i < numOfGenerations; i++) {
            
            //get population
            population = getPopulation(population);
            
            //handle crossover for population size
            population = handleCrossovers(population, nQeens);

            //loop until size to get fitness accuracy
            for (int j = 0; j < populationSize; j++) {

                if (getFitness(population[j]) == Fitness)
                    return population[j];
                
                //get proper mutation
                population[j] = mutate(population[j], mutation);
                if (getFitness(population[j]) == Fitness)
                    return population[j];
            }
        }
        return null;
    }
    private int[][] handleCrossovers(int[][] population, int n) {
        //loop until current population is reached
        for (int i = 0; i < population.length; i += 2) {
            
            //crossover Posibily formula
            int crossoverPos = (int) (Math.random() * n);

            //swap to get proper form of array[i][j]
            for (int j = 0; j < crossoverPos; j++) {
                int tmp = population[i][j];
                population[i][j] = population[i+1][j];
                population[i+1][j] = tmp;
            }
        }
        return population;
    }

    private int[][] getPopulation(int[][] population) {
        
        //Collections.sort(population, Comparator.comparingInt(this::getFitness));
        return population;
    }  
    //if small random probability then mutate, return result
    private int[] mutate(int[] res, double mutationProbability) {
        if (satisfyProb(mutationProbability))
            res[(int)(Math.random()*res.length)] = (int)(Math.random()*res.length);

        return res;
    }
    
    private boolean satisfyProb(double prob) {
        return prob >= Math.random();
    }

    private int getFitness(int[] res) {
        return getMaxFitness(res.length) - Solve.getHeuristicCost(res);
    }

    private int getMaxFitness(int n) {
        return n*(n-1)/2;
    }

    private int[] genIndividual(int n) {
        return Solve.setState(n);
    }

    private int[][] generatePopulation(int n, int populationSize) {
        int[][] population = new int[populationSize][];
        for (int i = 0; i < populationSize; i++)
            population[i] = genIndividual(n);

        return population;
    }
}
  class Solve {
    // Generate state that all queens have row # 0
    public static int[] oneDqueens(int n) {

        return new int[n];
    }
    // Randomizes state
    public static int[] randomState(int[] res) {

        for (int i = 0; i < res.length; i++)
            res[i] = (int) (Math.random() * res.length);

        return res;
    }
    // generates random initial state
    public static int[] setState(int n) {
        return randomState(oneDqueens(n));
    }
    // Returns heuristic cost
    public static int getHeuristicCost(int[] res) {
        int h = 0;
        // increment cost if two queens are in same row or in same diagonal.
        for (int i = 0; i < res.length; i++)
            for (int j = i + 1; j < res.length; j++)
                if (res[i] == res[j] || Math.abs(res[i] - res[j]) == j - i)
                    h += 1;
        return h;
    }
  }

class SAnnealing {

    public int[] solve(int boardSize, int maxNumOfIterations, double temperature, double coolingFactor) {
        int[] result = Solve.setState(boardSize);

        int costToBeat = Solve.getHeuristicCost(result);

        // terminate when it reaches max num of iterations or problem is solved.
        for (int x = 0; x < maxNumOfIterations && costToBeat > 0; x++) {
            result = move(result, costToBeat, temperature);
            costToBeat = Solve.getHeuristicCost(result);
            temperature = Math.max(temperature * coolingFactor, 0.01);
        }

        return costToBeat == 0 ? result : null; // return solution if solved
    }

    private int[] move(int res[], int costToBeat, double temp) {
        int n = res.length;

        while (true) {
            int i = (int) (Math.random() * n);
            int j = (int) (Math.random() * n);
            int tmpRow = res[i];
            res[i] = j;

            int cost = Solve.getHeuristicCost(res);
            if (cost < costToBeat)
                return res;

            int dE = costToBeat - cost;
            double acceptProb = Math.min(1, Math.exp(dE / temp));

            if (Math.random() < acceptProb)
                return res;

            res[i] = tmpRow;
        }
    }

}

public class NQueenProblem {
    public interface Solve {
public int[] solve();
} 
public static void printSolution(int[] res,int n, int board[][]){
if(res==null){
     System.out.print("No solution found ");
     }else 
{
      for (int i=0;i<n;i++){

          for (int j=0;j<n;j++){
              if(j!=res[i]){
              board[i][j]=0;
              }else
              board[i][j]=1;
          }
      }
 for (int i = 0; i < n; i++) { 
            for (int j = 0; j < n; j++) 
                System.out.print(" " + board[i][j] 
                                 + " "); 
            System.out.println(); 
        } 
}
}
public static void startGame(){
		System.out.println("Select an option to start N Queen Proble ");
		System.out.println("1. Simulated Annealing");
		System.out.println("2. GeneticAlgorithm");
		System.out.println("3. EXIT");
		System.out.print("\n->");
	};

    public static void main(String[] args) {
        
         int var1=0,var=0;
         SAnnealing simulatedAnnealing = new SAnnealing();
         GeneticAlgorithm geneticAlgorithm = new GeneticAlgorithm();
         int N=25;
         int saIterations=500;
         double saTemperature=120;
         double saCoolingFactor=.95;
         int populationSize=10;
         double mutationProbability=0.5;
         int numOfGenerations=50;
         int board[][]=new int[N][N];
        startGame();
		Scanner sc = new Scanner(System.in);
		var1=sc.nextInt(3);
		if (var1==1) {
                    long startTime = System.nanoTime();
                    int[] res = simulatedAnnealing.solve(N, saIterations, saTemperature, saCoolingFactor);
                    printSolution(res,N,board);
                    long endTime = System.nanoTime();
                    long timeElapsed = endTime - startTime;
                    System.out.println(" ");
                    System.out.println("Execution time in milliseconds : " + 
								timeElapsed / 1000000);
                    System.out.println(" ");
                    }
                else if(var1==2) {
                    long startTime = System.nanoTime();
                //int[] resGA= GeneticAlgorithm.solve(N,populationSize,mutationProbability,numOfGenerations);
                    //printSolution(resGA,N,board);
                    long endTime = System.nanoTime();
                    long timeElapsed = endTime - startTime;
                    System.out.println(" ");
                    System.out.println("Execution time in milliseconds : " + 
								timeElapsed / 1000000);
                    System.out.println(" ");
                }
                else if(var1==3) {
			System.exit(0);
		}
		sc.close();     
    }
    
}
