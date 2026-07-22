# CoreForge — Doğrulama Raporu (Simülasyon Kod Doğrulama)

Teknofest NDK şartnamesi gereği, yazılım **uluslararası kabul görmüş
benchmark verileriyle** test edilmiştir. Bu rapordaki HER sayı, bu kod
sistemiyle üretilmiş gerçek koşu sonucudur ve `python3 verify.py`
komutuyla herkes tarafından yeniden üretilebilir (24 otomatik kontrol).

## 1 · Doğrulama metodolojisi

İki kavram ayrıştırılır:
- **Verification (doğru çözüyor muyuz?):** analitik çözümlü problemler +
  mesh yakınsaması + değişmezlik testleri.
- **Validation (doğru modeli mi çözüyoruz?):** uluslararası benchmark
  referanslarıyla karşılaştırma; model sınırının kendisi de ölçülür
  (C5G7 vakası).

## 2 · Analitik doğrulamalar (kesin referans)

| Vaka | Model | Sonuç | Analitik | Fark |
|---|---|---|---|---|
| Homojen sonsuz ortam k∞ | 2-D | 1.0857142 | 1.0857143 | −0.0 pcm |
| Çıplak kare, 1 grup, B²=2(π/L)² | 2-D | 1.2610188 | 1.2610203 | −0.1 pcm |
| Çıplak küp, 1 grup, B²=3(π/L)² | 3-D | 1.2490752 | 1.2490256 | +3.2 pcm |
| Nokta kinetiği, +100 pcm periyodu | transient | 54.85 s | inhour 54.92 s | %0.12 |
| Nokta kinetiği, −100 pcm periyodu | transient | −130.0 s | inhour −130.0 s | %0.02 |
| T-H enerji kapanışı (tek kanal) | ısıl-hidrolik | 4.3200 MW | 4.32 MW | tam |
| Yakıt ΔT el-hesabı (20 kW/m) | ısıl-hidrolik | 531 K | q'/4πk = 530.5 K | ✓ |
| İyot çukuru zamanı/derinliği | xenon | 7.6 h / ×1.55 | PWR bandı 6–13 h | ✓ |

## 3 · Uluslararası benchmark doğrulamaları

### IAEA-2D PWR (ANL-7416, problem 11-A2) — referans k=1.02959
| h [cm] | k_eff | Fark |
|---|---|---|
| 2 | 1.0295069 | −7.8 pcm |
| 1 | 1.0295843 | **−0.5 pcm** |
| 0.5 | 1.0296051 | +1.4 pcm |
| h→0 (Richardson) | 1.0295975 | **+0.7 pcm, order 1.90** |

Monoton h² yakınsaması (v8.4 düzgün-yakınsamış çözücü ile order 1.90);
geometri resmî poligon tanımından 10 cm blok kafesine birebir
aktarılmıştır. (v8.3'teki "−1.1 pcm" değeri kısmen yanlış-yakınsama
artıfaktıydı; bkz. §çözücü notu.)

### IAEA-3D PWR (ANL-7416, problem 11) — referans k=1.02903
Tam x-y-z: 380 cm, eksenel reflektörler, 4 tam çubuk + üstten 80 cm
batmış 5. çubuk (kutu koordinatları resmî tanımdan doğrulanmıştır).

| div / divz | k_eff | Fark |
|---|---|---|
| 2 / 2 | 1.028666 | −34.4 pcm |
| 4 / 4 | 1.028896 | **−12.7 pcm** (referansa doğru monoton) |

(v8.3'teki "−6.7 pcm" değeri yanlış-yakınsama artıfaktıydı; v8.4 düzgün-
yakınsamış çözücüyle gerçek değer −12.7 pcm ve mesh inceldikçe referansa
monoton yaklaşır.)

Eksenel güç profili beklenen üst-bastırılmış asimetriyi verir
(F_z=1.56). **Çubuk-5 batırma taraması**: toplam değer ~−980 pcm,
klasik S-eğrisi (ilk 85 cm'de −85 pcm, orta korda −476 pcm).

### C5G7 MOX (OECD/NEA, NEA/NSC/DOC(2003)16) — transport referansı 1.18655
Pin-hücre homojenize **difüzyon demosu** (kesitler doğrulanmış S_N
transport çözücümüzden birebir):

| h [cm] | k_eff | Fark |
|---|---|---|
| 1.26 | 1.184796 | −125 pcm |
| 0.63 | 1.186392 | −11 pcm |
| h→0 (Richardson) | 1.186984 | **+31 pcm = saf model hatası** |

Özdeğer uyumu kısmen hata iptalidir; pin güçleri sapar (maks. pin 2.42
vs S_N 2.35) — difüzyonun pin-ölçeği sınırının NİCEL kanıtı.

## 4 · Değişmezlik ve tutarlılık testleri

| Test | Sonuç |
|---|---|
| 2-D çözüm ≡ tek-zon 3-D (yansıtıcı alt/üst) | 1.0292601 vs 1.0292600 → **0.01 pcm** |
| Taze yakıt: cell_from_N ≡ pincell_xs | bit-düzeyi özdeş |
| Nötron dengesi ↔ analitik kaçak (çıplak kor) | %3.60 vs %3.64 |
| Burnup enerji muhasebesi | yanma haritası ort. = hedef (tam) |
| Mutlak akı kalibrasyonu (enerji kapanışı) | 50.000 MW ↔ 50 MW hedef (1-grup testi) |
| QA parmakizi (bütünlük kilidi) | mesh inceltme referansı KORUR; pitch/kesit/BC değişimi otomatik GEÇERSİZ kılar |

## 5 · Fizik trend doğrulamaları (yakıt tasarımcısı & çevrim)

| Büyüklük | Sonuç | Literatür bandı |
|---|---|---|
| k∞(zenginlik) 2.1/3.1/4.5 w/o | 1.217/1.320/1.408 (monoton) | ✓ |
| Boron değeri | −6.1 pcm/ppm | −5…−11 |
| Bağımsız kaynak çaprazı (IAEA yakıtı) | en kötü sabit farkı %17.5 (Σ₁₂) | <%25 |
| Ters çözücü: IAEA yakıt-1 | e=2.886 w/o + 259 ppm (νΣf₂/Σa₂ tam) | — |
| Reaktivite-sınırlı yanma (%3.1) | 39.7 MWd/kgU | 28–45 |
| Xe+Sm denge değeri (kor) | −3 567 pcm | −2 500…−4 000 |
| Kritik boron (taze, zehirsiz) | 2 236 ppm → k=1.000000 | taze kor sınıfı |
| Boron letdown | 1 767→0 ppm, her adım k=1.0000 | monoton |
| Otomatik çevrim sonu | 30.0 MWd/kgU = 789 EFPD | — |
| Çevrimde güç düzleşmesi | F_xy 2.60 → 1.32 | klasik |
| REA (rod ejection, trip'li) | scram güç %4.2'ye kapatır | klasik |
| ATWS (scram başarısız) | Doppler+MTC %167'de dengede tutar | scram değeri görünür |

## 6 · Bağımsız kod çaprazı

Aynı ekipçe geliştirilen **ayrı S_N transport çözücüsü**
(github.com/EmreSakarya/c5g7-2d-transport-benchmark, MCNP'ye −182 pcm)
ile C5G7 üzerinde kod-koda karşılaştırma yapılmıştır; ayrıca IAEA-3D,
ekibin bağımsız 3-D difüzyon kodunun sonucuyla aynı mertebede
(v8.4 düzgün-yakınsamış değer −12.7 pcm) doğrulanmıştır.

### Çözücü notu (v8.4 — yanlış-yakınsama düzeltmesi)

v8.3'te sabit iç-SOR süpürme sayısı ince mesh'te flux'u yakınsatmıyor,
dış Δk-testi yavaş-sürünmede yanlış `converged` verebiliyordu (mesh
inceldikçe k referanstan uzaklaşırdı). v8.4 iç çözümü **kaynak-normalize
rezidu** ile adaptif yapar; mesh yakınsaması artık monoton order-2 ve
`verify.py`'deki `convergence_check` bunu kalıcı olarak kilitler.
Varsayılan-mesh ve C5G7 değerleri değişmedi; yalnız ince-mesh/Richardson
değerleri düzgün-yakınsamış hâle getirildi.

## 7 · Yeniden üretilebilirlik

```bash
python3 verify.py              # 26 kontrol (analitik+benchmark+trend)
python3 verify.py --fine       # ince-mesh satırları dâhil
python3 verify.py --no-engine  # Fortran'sız saf-Python fizik alt kümesi
```
Her arayüz çözümü, girdi dosyası ve CSV'leriyle birlikte indirilebilir;
📄 HTML raporu tüm konfigürasyon ve sonuçları tek dosyada arşivler.
