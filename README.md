# Trabalho de Cálculo Numérico (Gêmeo Digital)

Projeto da disciplina SME0602 - Cálculo Numérico (ICMC-USP), com a construção incremental de um Gêmeo Digital de um dispositivo microfluídico, através da aplicação de métodos numéricos a seis subsistemas físicos acoplados.

## Instalação

```bash
pip install -r requirements.txt
```

## Como executar

Cada capítulo é independente e roda a partir da sua própria pasta:

```bash
python NomeDaPasta/main*.py
```

## Estrutura

| Capítulo | Pasta | Arquivos |
|---|---|---|
| 1. Rede Hidráulica | `RedeHidraulica/` | `main.py`, `functions.py` |
| 2. Placa Térmica | `PlacaTermica/` | `mainT.py`, `functionsT.py` |
| 3. Membrana Elástica | `MembranaElastica/` | `mainM.py`, `functionsM.py` |
| 4. Acoplamento Hidráulico-Térmico | `AcoplamentoHidraulicoTermico/` | `mainHT.py`, `functionsHT.py` |
| 5. Acoplamento Hidráulico-Mecânico | `AcoplamentoHidraulicoMecanico/` | `mainHM.py`, `functionsHM.py` |
| 6. Gêmeo Digital | `GemeoDigital/` | `mainGD.py`, `functionsGD.py` |

As pastas também possuem outros arquivos extras além do código, como as apresentações, algumas imagens dos plots e algumas tabelas dos resultados

## Relatório

O relatório final, com a fundamentação teórica e a análise de todos os resultados, está disponível no arquivo `RelatorioGrupo1.pdf`.
