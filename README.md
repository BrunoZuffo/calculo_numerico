# calculo_numerico
# 🐍 Setup do Ambiente Python (venv + requirements)

Este projeto utiliza **ambiente virtual (venv)** para isolar dependências.  
Cada integrante deve criar sua própria venv localmente.

⚠️ A pasta `.venv/` **não deve ser enviada para o Git**.

---

## 📦 1️⃣ Clonar o projeto

```bash
git clone <URL_DO_REPOSITORIO>
cd <NOME_DO_PROJETO>
🐧 2️⃣ Setup no Linux / Ubuntu / Mac
1. Criar a venv
python3 -m venv .venv
2. Ativar a venv
source .venv/bin/activate

Você saberá que está ativa quando aparecer (.venv) no terminal.

3. Instalar dependências
pip install -r requirements.txt
🪟 3️⃣ Setup no Windows (PowerShell)
1. Criar a venv
python -m venv .venv
2. Ativar a venv
.venv\Scripts\Activate

Se der erro de permissão, rode antes:

Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

Depois tente ativar novamente.

3. Instalar dependências
pip install -r requirements.txt
🔄 Sempre que abrir o projeto

Você precisa ativar a venv novamente.

Linux:

source .venv/bin/activate

Windows:

.venv\Scripts\Activate
➕ Quando adicionar uma nova biblioteca

Se você instalar algo novo:

pip install nome-da-biblioteca

Atualize o requirements.txt:

pip freeze > requirements.txt

Depois:

git add requirements.txt
git commit -m "Update requirements"
git push
🧹 Arquivos que NÃO devem ir para o Git

Adicione ao .gitignore:

.venv/
__pycache__/
*.pyc

Nunca envie sua venv para o repositório.


---

Se quiser, posso agora:

- Deixar isso mais profissional estilo projeto open-source  
- Ou mais enxuto para projeto acadêmico  
- Ou adicionar seção de execução do projeto (`python main.py`, etc.)