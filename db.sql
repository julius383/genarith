CREATE TABLE IF NOT EXISTS experiment_stats (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  target REAL UNIQUE NOT NULL,
  population_size INTEGER NOT NULL,
  max_generations INTEGER NOT NULL,
  generations_run INTEGER NOT NULL,
  mutations INTEGER NOT NULL,
  crossovers INTEGER NOT NULL,
  solution TEXT,
  population JSON
);
