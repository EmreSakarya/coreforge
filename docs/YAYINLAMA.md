# Yayınlama Kılavuzu — GitHub + Streamlit Cloud Canlı Demo

Bu doküman CoreForge'u (1) GitHub'da yayınlamak ve (2) Streamlit
Community Cloud'da ücretsiz canlı demo olarak açmak için gereken TÜM
adımları içerir. Komutlar WSL/Ubuntu terminali içindir.

---

## 0. Güvenlik kuralları (önce oku)

* **Token'ı asla** bir sohbete, dosyaya, ekran görüntüsüne veya commit'e
  yapıştırmayın. Token yalnızca `git remote` URL'sinde ya da git'in
  kimlik doğrulama penceresinde kullanılır.
* Yanlışlıkla paylaşılan bir token'ı hemen iptal edin:
  GitHub → Settings → Developer settings → Personal access tokens →
  ilgili token → **Revoke**.
* Depoya `runs/`, derlenmiş motor (`solver/coreforge`) ve `*.mod`
  dosyaları girmez — `.gitignore` bunu zaten sağlıyor; yine de ilk
  push'tan önce `git status` çıktısını kontrol edin.

## 1. GitHub deposunu oluştur

Tarayıcıda: github.com → **New repository**

* Repository name: `coreforge`
* Description: `2D/3D multigroup neutron diffusion reactor analysis code — Fortran engine + Streamlit UI (IAEA-2D/3D, C5G7 validated)`
* **Public** (Streamlit Community Cloud ücretsiz katmanı public depo ister)
* README/License/gitignore eklemeyin (hepsi pakette hazır) → **Create**

## 2. Yerel depoyu hazırla ve gönder

Proje klasöründe (`coreforge/`):

```bash
cd coreforge

# kimlik — depo katkıcısı yalnızca siz olacaksınız
git config --global user.name  "Emre Sakarya"
git config --global user.email "gemrek123@gmail.com"

git init -b main
git add .
git status          # runs/, solver/coreforge, __pycache__ LİSTEDE OLMAMALI
git commit -m "CoreForge v8.4: 2D/3D multigroup diffusion reactor analysis code"

git remote add origin https://github.com/EmreSakarya/coreforge.git
git push -u origin main
```

**Kimlik doğrulama (token güvenliği).** `git push` sizi token/parola
sorar. Token'ı `git remote` URL'sine gömerseniz **düz metin olarak
`.git/config`'te kalıcılaşır** — bunun yerine güvenli yolu kullanın:

* En temizi: `gh auth login` (GitHub CLI) veya bir kimlik yardımcısı:
  ```bash
  git config --global credential.helper store   # ya da: cache
  ```
  İlk push'ta kullanıcı adı + token isteyince token'ı **parola alanına**
  yapıştırın; helper güvenli biçimde saklar.
* Token'ı komut satırına, bir dosyaya, ekran görüntüsüne veya bir sohbete
  **asla** yapıştırmayın. Yanlışlıkla ifşa olursa hemen GitHub →
  Settings → Developer settings'ten **Revoke** edin.
* URL'ye gömmek zorunda kalırsanız (`https://<TOKEN>@github.com/...`),
  push'tan sonra `git remote set-url origin
  https://github.com/EmreSakarya/coreforge.git` ile token'ı config'ten
  temizleyin.

Push sonrası kontrol listesi:

- [ ] github.com/EmreSakarya/coreforge açılıyor, README ekran
      görüntüleriyle birlikte görünüyor
- [ ] Sağ panelde "MIT License" rozeti var
- [ ] Contributors listesinde yalnızca **Emre Sakarya** var
- [ ] `runs/`, `solver/coreforge` (binary), `__pycache__` depoda YOK

Depoya "topics" ekleyin (About → ⚙): `nuclear-engineering`,
`reactor-physics`, `neutron-diffusion`, `fortran`, `streamlit`,
`benchmark`, `smr`. Keşfedilebilirliği ciddi artırır.

## 3. Streamlit Community Cloud canlı demo

Paket bulut için hazır: `requirements.txt` Python bağımlılıklarını,
`packages.txt` sistem paketi `gfortran`'ı kurdurur; uygulama ilk
açılışta Fortran motorunu **kendisi derler** (`runner.build_engine`),
derleme günlüğü `solver/build.log`'a düşer.

1. https://share.streamlit.io → **Sign in with GitHub**
2. **Create app** → *Deploy a public app from GitHub*
3. Repository: `EmreSakarya/coreforge` · Branch: `main` ·
   Main file path: `app.py`
4. App URL alt alanını seçin (ör. `coreforge`) → **Deploy**

İlk açılış 2–4 dakika sürer (paket kurulumu + motor derlemesi).
Uygulama açılınca **Benchmarks** sekmesinden IAEA-2D'yi (varsayılan
mesh, h=2 cm) çalıştırıp k≈1.02951 (referansa −8 pcm) sonucunu görerek
bulut kurulumunu doğrulayın; mesh'i inceltirseniz (kenar çubuğu → div)
sonuç −1 pcm'e yakınsar.

Notlar:

* Ücretsiz katman ~1 GB RAM verir: varsayılan çözünürlükler rahat
  çalışır; çok ince mesh (div≥8 + divz≥4 3B) bulutta yavaş olabilir —
  bu normaldir, README'de "canlı demo hafif ayarlar içindir" notu var.
* Uygulama bir süre ziyaret edilmeyince uyur; ilk ziyaretçi geldiğinde
  kendiliğinden uyanır (motor tekrar derlenmez, kalıcı disk aynı kaldığı
  sürece).
* Canlı demo linkini README'nin en üstüne ekleyin:
  `**Canlı demo:** https://<secilen-ad>.streamlit.app`

## 4. Sonraki güncellemeler

```bash
git add -A
git commit -m "kısa ve açıklayıcı mesaj"
git push
```

`main`'e her push, Streamlit Cloud uygulamasını otomatik yeniden dağıtır.

## 5. Sürümleme önerisi

Kararlı her aşamada etiket atın:

```bash
git tag -a v8.4 -m "CoreForge v8.4 — publication release"
git push origin v8.4
```

GitHub → Releases → **Draft a new release** ile v8.4 etiketinden bir
sürüm yayınlayıp `coreforge.zip`'i eklerseniz, kullanıcılar tek tıkla
indirilebilir arşive kavuşur (Teknofest başvurusunda da temiz bir
"prototip sürümü" referansı olur).
