#include <random>
#include <vector>

struct HardDisks
{
  enum class Method
  {
    Dynamic,
    Static
  };

  size_t numberOfInitCycles;
  size_t numberOfProdCycles;
  size_t numberOfParticles;
  double maxDisplacement;
  size_t sampleFrequency;
  double boxSize;

  Method method;

  std::vector<std::pair<double, double>> positions;

  size_t numberOfAttemptedMoves;
  size_t numberOfAcceptedMoves;

  std::random_device rd;
  std::mt19937 gen(rd());

  HardDisks(size_t numberOfInitCycles, size_t numberOfProdCycles, size_t numberOfInitParticles, double maxDisplacement,
            size_t sampleFrequency, double boxSize, Method method);

  void initialize();
  void run();
  void getMSD();
}