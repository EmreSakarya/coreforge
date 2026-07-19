# CoreForge — Yayınlama Kılavuzu (GitHub + Streamlit Cloud)

Bu kılavuz, CoreForge'u herkese açık bir GitHub deposu olarak yayınlamak ve
Streamlit Community Cloud'da canlı demoya almak için gereken adımları içerir.
Depo halihazırda yayına hazırdır: `LICENSE` (MIT), `requirements.txt`,
`packages.txt` (bulutta `gfortran` kurar), motor otomatik-derleme (`runner.py`
ilk çağrıda derler) ve `.streamlit/config.toml` mevcuttur.

> **Güvenlik kuralı:** GitHub kişisel erişim token'ınızı (PAT) hiçbir sohbete,
> dosyaya, commit'e veya ekran görüntüsüne **yapıştırmayın**. Aşağıdaki iki
> yöntem de token'ı elle taşımanızı gerektirmez: `gh auth login` tarayıcı/
> cihaz akışı kullanır, SSH ise anahtar tabanlıdır.

---

## 0 · Ön koşullar

- Yerelde ilk commit yapılmış olmalı (kimlik: `Emre Sakarya
  <gemrek123@gmail.com>`, dal: `main`). Kontrol:
  ```bash
  git log -1 --format='%h %an <%ae>'   # Emre Sakarya görünmeli
  git branch --show-current            # main
  ```
- GitHub hesabı: kullanıcı adı sizin (repo adı önerisi: `coreforge`, **public**).

---

## 1 · Push — Yöntem A: GitHub CLI (`gh`) — önerilen, token yapıştırma yok

`gh auth login`, tarayıcı üzerinden tek-kullanımlık bir cihaz kodu ile
kimlik doğrular; token'ı siz görmezsiniz, hiçbir yere yapıştırmazsınız.

**gh kurulumu (WSL2 / Ubuntu):**
```bash
sudo apt update && sudo apt install -y gh
# veya resmi depo: https://github.com/cli/cli/blob/trunk/docs/install_linux.md
```

**Giriş (interaktif — bu komutu kendiniz terminalde çalıştırın):**
```bash
gh auth login
#   ? What account do you want to log into?  GitHub.com
#   ? Preferred protocol for Git operations?  HTTPS
#   ? Authenticate Git with your GitHub credentials?  Yes
#   ? How would you like to authenticate?  Login with a web browser
#   -> ekranda bir tek-kullanımlık kod görünür; Enter'a basınca tarayıcı
#      açılır, kodu girip onaylarsınız. Token otomatik saklanır.
```

**Depoyu oluştur ve push et (tek komut):**
```bash
gh repo create coreforge --public --source=. --remote=origin --push
```
Bu komut GitHub'da `coreforge` deposunu açar, `origin` remote'unu ekler ve
`main` dalını push eder. Bitti.

---

## 2 · Push — Yöntem B: SSH anahtarı — token tamamen yok

Token hiç kullanmak istemiyorsanız SSH en temiz yoldur.

```bash
# 1) Anahtar üret (zaten varsa atla: ls ~/.ssh/id_ed25519.pub)
ssh-keygen -t ed25519 -C "gemrek123@gmail.com"    # Enter'larla varsayılanları kabul et

# 2) Public anahtarı kopyala ve GitHub'a ekle:
cat ~/.ssh/id_ed25519.pub
#   -> GitHub > Settings > SSH and GPG keys > New SSH key > yapıştır
#      (bu anahtardır, TOKEN DEĞİL — public anahtar paylaşmak güvenlidir)

# 3) Bağlantıyı test et
ssh -T git@github.com    # "Hi <kullanıcı>! You've successfully authenticated" görmelisin
```

Ardından GitHub arayüzünden boş, **README'siz** bir public `coreforge` deposu
açın (ekstra dosya eklemeyin ki push çakışmasın), sonra:
```bash
git remote add origin git@github.com:<KULLANICI_ADI>/coreforge.git
git push -u origin main
```

---

## 3 · Push sonrası doğrulama

```bash
git remote -v                 # origin görünmeli
git log origin/main -1        # commit uzakta olmalı
```
GitHub'da depoyu açın; `README.md` otomatik render olmalı, `docs/` içindeki
ekran görüntüleri görünmeli. Motor binary'si, `runs/` ve `__pycache__`
`.gitignore` sayesinde **yüklenmemiş** olmalı (bu kasıtlıdır — bulut kendi
derler).

---

## 4 · Streamlit Community Cloud — canlı demo

Depo public olunca:

1. <https://share.streamlit.io> → GitHub ile giriş yapın.
2. **New app** → repo: `<KULLANICI_ADI>/coreforge`, branch: `main`,
   main file: `app.py`.
3. **Deploy**. İlk açılışta:
   - `packages.txt` içindeki `gfortran` sistem paketini kurar,
   - `requirements.txt` Python bağımlılıklarını kurar,
   - motor ilk çalıştırmada `runner.py` tarafından bulutta **otomatik
     derlenir** (log: `solver/build.log`). ifx bulutta yoktur; gfortran
     yedeği devreye girer.
4. Birkaç dakika sonra `https://<app-adı>.streamlit.app` adresinde canlı olur.

**İpuçları**
- İlk soğuk açılış motoru derlediği için birkaç saniye uzun sürer; sonraki
  koşular hızlıdır.
- App güncellemek için sadece `git push` yeterlidir; Streamlit otomatik
  yeniden dağıtır.

---

## 5 · Yayın öncesi son kontrol listesi

- [ ] `python3 verify.py --fine` → tüm kontroller PASS (sürüm-lockstep dahil).
- [ ] `git log -1 --format='%an <%ae>'` → `Emre Sakarya <gemrek123@gmail.com>`.
- [ ] `git ls-files` çıktısında binary / `runs/` / `__pycache__` **yok**.
- [ ] `README.md` ve `LICENSE` mevcut; sürüm dizeleri (motor/report/app) 8.3.
- [ ] Repo **public**, ad `coreforge`.
