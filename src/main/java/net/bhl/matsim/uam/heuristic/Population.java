package net.bhl.matsim.uam.heuristic;

import net.bhl.matsim.uam.analysis.traveltimes.utils.TripItem;
import net.bhl.matsim.uam.optimization.Vertiport;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class Population {

        private Chromosome[] population;
        public List<Vertiport> vertiportCandidates;
        public List<TripItem> deserializedTripItems;

        public Population(int populationSize,List<Vertiport> vertiportCandidates, List<TripItem> deserializedTripItems) {
            this.population = new Chromosome[populationSize];
            this.vertiportCandidates=vertiportCandidates;
            this.deserializedTripItems=deserializedTripItems;
        }

        public void initializePopulation() {
            for (int i = 0; i < this.population.length; i++) {
                this.population[i] = new Chromosome(generateRandomGenes(200));
            }
        }

        public void calculateFitnessOfAllChromosomes() {
            for (Chromosome chromosome : this.population) {
                chromosome.calculateFitness(this.vertiportCandidates,this.deserializedTripItems);
            }
        }

        public List<Integer> generateRandomGenes(int length) {
          // generate a Integer List with a length of 200, 74 of them are 1, the rest are 0
            List<Integer> genes = new ArrayList<>();
            for (int i = 0; i < length; i++) {
                if (i < 74) {
                    genes.add(1);
                } else {
                    genes.add(0);
                }
            }
            Collections.shuffle(genes); // TODO: check if this is a good way to generate random genes, maybe use random class to keep it repeatable
            return genes;
        }

        public void select() {
            // 锦标赛选择
            List<Chromosome> selected = new ArrayList<>();
            int tournamentSize = 5;

            for (int i = 0; i < this.population.length; i++) {
                List<Chromosome> tournament = new ArrayList<>();
                for (int j = 0; j < tournamentSize; j++) {
                    tournament.add(this.population[new Random().nextInt(this.population.length)]);
                }

                // sort the tournament by fitness in descending order
                tournament.sort((o1, o2) -> {
                    if (o1.fitness > o2.fitness) {
                        return -1;
                    } else if (o1.fitness < o2.fitness) {
                        return 1;
                    } else {
                        return 0;
                    }
                });
                // add the fittest chromosome in the tournament to the selected list
                selected.add(tournament.get(0));
            }

            for (int i = 0; i < this.population.length; i++) {
                this.population[i] = selected.get(i);
            }
        }

        public void crossover() {
            // 单点交叉
            for (int i = 0; i < this.population.length - 1; i += 2) {
                int crossoverPoint = new Random().nextInt(200);
                List<Integer> parent1 = new ArrayList<>(this.population[i].genes);
                List<Integer> parent2 = new ArrayList<>(this.population[i + 1].genes);
                List<Integer> child1 = new ArrayList<>(parent1);
                List<Integer> child2 = new ArrayList<>(parent2);
                for (int j = crossoverPoint; j < 200; j++) {
                    int temp = parent1.get(j);
                    child1.set(j, parent2.get(j));
                    child2.set(j, temp);
                }
                // check the validity of the child: the number of 1 should be 74, if not, change some code of the not changed part
                while (child1.stream().mapToInt(Integer::intValue).sum() < 74) {
                  Random random = new Random();
                    int index = random.nextInt(crossoverPoint);
                    if (child1.get(index) == 0) {
                        child1.set(index, 1);
                    }
                }
                while (child1.stream().mapToInt(Integer::intValue).sum() > 74) {
                    Random random = new Random();
                    int index = random.nextInt(crossoverPoint);
                    if (child1.get(index) == 1) {
                        child1.set(index, 0);
                    }
                }

                while (child2.stream().mapToInt(Integer::intValue).sum() < 74) {
                    Random random = new Random();
                    int index = random.nextInt(crossoverPoint);
                    if (child2.get(index) == 0) {
                        child2.set(index, 1);
                    }
                }
                while (child2.stream().mapToInt(Integer::intValue).sum() > 74) {
                    Random random = new Random();
                    int index = random.nextInt(crossoverPoint);
                    if (child2.get(index) == 1) {
                        child2.set(index, 0);
                    }
                }

                this.population[i]=new Chromosome(child1);
                this.population[i+1]=new Chromosome(child2);
            }
        }

        public void mutate() {
            // randomly change one 0 to 1
            for (Chromosome chromosome : this.population) {
                if (new Random().nextDouble() < 0.05) {
                    int index1 = new Random().nextInt(200);
                    while (chromosome.genes.get(index1) == 1) {
                        index1 = new Random().nextInt(200);
                    }
                    chromosome.genes.set(index1, 1);

                    int index2 = new Random().nextInt(200);
                    while (chromosome.genes.get(index2) == 0) {
                        index2 = new Random().nextInt(200);
                    }
                    chromosome.genes.set(index2, 0);
                }
            }

        }

        public Chromosome getbestChromosome() {
            Chromosome fittest = this.population[0];
            for (Chromosome chromosome : population) {
                if (chromosome.fitness > fittest.fitness) {
                    fittest = chromosome;
                }
            }
            return fittest;
        }

    }

