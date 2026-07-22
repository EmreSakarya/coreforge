# Bağımsız Doğrulama Kılavuzu (Independent Verification Pack)

Bu doküman, CoreForge'un doğruluğunu ve kalitesini **bu projeyi geliştiren
araçtan bağımsız** bir yapay zekâya (veya bir insana) test ettirmek için
hazır malzeme içerir. Üç tur önerilir; her turun kopyala-yapıştır
kullanılabilecek istemi (prompt) aşağıdadır.

Neden bağımsız doğrulama? Bir sistemi geliştiren araç, kendi kör
noktalarını paylaşır. Bağımsız bir gözün aynı kanıtları yeniden üretmesi
(veya üretememesi) en güçlü kalite sinyalidir. Bu, nükleer yazılım
pratiğindeki bağımsız gözden geçirme (independent review) ilkesinin
küçültülmüş halidir.

---

## Hazırlık

Bağımsız doğrulayıcıya verilecek malzeme:

1. **Kod:** `coreforge.zip` dosyasının tamamı (veya GitHub deposu linki).
2. **İddialar:** `README.md` içindeki doğrulama tabloları ve
   `docs/DOGRULAMA_RAPORU.md`.
3. Aşağıdaki istemlerden uygun olanı.

> Not: Dosya yükleme destekleyen ve kod çalıştırabilen bir yapay zekâ
> seçin (ör. kod yorumlayıcısı/sandbox özelliği olan bir sohbet modeli).
> Sandbox'ta Fortran derleyicisi yoksa `--no-engine` yolu kullanılır
> (Tur A'da açıklanıyor).

---

## Tur A — Çalıştırmalı doğrulama (en güçlü kanıt)

Amaç: 24 otomatik kontrolün bağımsız bir makinede aynı sonuçları
üretmesi. Kopyalayıp yapıştırın:

```text
Ekteki coreforge.zip bir 2B/3B çok gruplu nötron difüzyon çözücüsü
(Fortran motoru + Python sürücüsü). Benim yazdığım bir kod DEĞİL;
tarafsız biçimde doğrulamanı istiyorum.

1) Zip'i aç. Ortamında Fortran derleyicisi varsa motoru derle:
     gfortran -O3 -fopenmp solver/coreforge.f90 -o solver/coreforge
   (Python tarafı zaten eksik motoru gfortran ile kendisi de derlemeyi
   dener.) numpy gereklidir: pip install numpy
2) python3 verify.py  çalıştır ve çıkan tabloyu OLDUĞU GİBİ raporla.
   Fortran derleyicin yoksa şunu çalıştır:
     python3 verify.py --no-engine
   Bu, motorsuz saf-Python fizik alt kümesidir (tasarımcı eğilimleri,
   tükenme 0-B, termal-hidrolik el hesapları, ksenon, kinetik).
3) Şu üç soruya kanıta dayalı cevap ver:
   a) PASS/FAIL tablosu koddaki iddialarla (README doğrulama tablosu)
      tutarlı mı?
   b) Kontrollerden herhangi biri "kendi kendini doğrulayan" (tautolojik)
      mi — yani test, test ettiği kodla aynı formülü mü kullanıyor?
      Hangileri gerçekten bağımsız referanslara (IAEA/OECD benchmark
      değerleri, el hesapları, analitik çözümler) dayanıyor?
   c) Testlerin KAPSAMADIĞI önemli bir fiziksel davranış görüyor musun?
```

Beklenen sonuç: `all checks passed` (motorlu yolda 26 kontrol,
motorsuz yolda 13 kontrol). Tek bir FAIL bile önemlidir — bana bildirin.

---

## Tur B — Düşmanca (adversarial) fizik ve kod eleştirisi

Amaç: hata bulmaya odaklı bağımsız bir göz. Model kapsamı bilinçli
yaklaşımları "hata" sanmasın diye kapsam bildirimi isteme dahildir.
Kopyalayıp yapıştırın:

```text
Ekteki coreforge.zip'teki reaktör fiziği kodunu DÜŞMANCA incele: görevin
onu övmek değil, kusur bulmak. Kapsam bildirimi: bu kasıtlı olarak bir
difüzyon kodudur (transport değil), 2 grup yarı-ampirik XS üretir
(ENDF işleme değil), kapalı-kanal tek-faz T-H yapar (alt-kanal kodu
değil). Bu kapsam SINIRLARINI hata olarak sayma; kapsam İÇİNDEKİ
hataları ara.

Modül modül incele ve her bulgu için (dosya:satır, ciddiyet, somut
hatalı senaryo) ver:
1) solver/coreforge.f90 — FDM ayrıklaştırma, harmonik ortalama D,
   Robin vakum sınır koşulu (J = 0.4692·φ), kırmızı-siyah SOR,
   güç iterasyonu. Şablon katsayılarında işaret/faktör hatası var mı?
2) xslib.py — 2200 m/s mikroskopikler doğru mu? (σf5=582.6 b,
   σa5=680.9 b, ν=2.432, σaB=759 b, Xe σa=2.65e6 b ...) Spektrum
   faktörü yaklaşımının bozulacağı durumlar hangileri?
3) burnup.py — Bateman RK4, denge Xe/Sm, bor letdown. Birim hataları
   (barn↔cm², J↔MeV, gün↔saniye) var mı?
4) kinetics.py — Keepin 6 grup (β=0.006502), inhour denklemi, matris
   üsteli. Negatif reaktivitede kök seçimi doğru mu?
5) thermal.py — Dittus-Boelter, iletim dirençleri, doyma marjı.
   Enerji korunumu gerçekten kapanıyor mu?
6) app.py/runner.py — durum yönetimi, fingerprint bütünlük katmanı
   atlatılabilir mi (referans geçersizken geçerli gösterilebilir mi)?

En ciddi 5 bulgunu ciddiyet sırasıyla listele. Bulgu yoksa "yok" deme;
neden bulamadığını (hangi kontrolün seni ikna ettiğini) açıkla.
```

Değerlendirme: MAJOR bulgu → düzeltilmeli ve verify.py'ye kalıcı kontrol
eklenmeli (bu projenin standart pratiği). MINOR/kapsam-dışı → yol
haritasına not.

---

## Tur C — Benchmark referanslarının literatür kontrolü

Amaç: kodun KENDİSİNİ değil, karşılaştırıldığı REFERANSLARI doğrulamak.
Web erişimi olan bir yapay zekâya kopyalayıp yapıştırın:

```text
Aşağıdaki referans değerlerin nükleer mühendislik literatüründeki
kaynaklarını bul ve doğrula; her biri için kaynak (rapor no/DOI) ver:

1) IAEA-2D PWR benchmark (ANL Benchmark Problem Book, problem 11-A2):
   k_eff referansı 1.02959 mu?
2) IAEA-3D PWR benchmark (problem 11): k_eff referansı 1.02903 mü?
   (FeenoX ve ADPRES/KOMODO dokümanlarında da alıntılanır.)
3) OECD/NEA C5G7 2-D MOX benchmark (NEA/NSC/DOC(2003)16): MCNP
   referansı k_eff = 1.18655 mi?
4) Keepin U-235 termal fisyon gecikmiş nötron verileri: toplam
   β = 0.0065 civarı mı; 6 grup λ değerleri hangi kaynakta?
5) U-235 2200 m/s kesitleri: σ_f ≈ 582.6 b, σ_a ≈ 680.9 b (ENDF/B ve
   Atlas of Neutron Resonances ile uyumlu mu)?
6) 15.5 MPa'da su doyma sıcaklığı ≈ 344.8 °C (IAPWS-IF97) doğru mu?

Ayrıca: bu benchmark'lar için literatürde raporlanmış TİPİK difüzyon
kodu sapmaları nedir? (Ör. C5G7'de difüzyonun yüzlerce pcm sapması
normal midir?) Ekteki README'nin doğrulama tablosundaki sapmalar bu
tipik aralıklarla tutarlı mı?
```

---

## Sonuçların yorumlanması

| Sonuç | Anlamı | Aksiyon |
|---|---|---|
| Tur A tabloları birebir aynı | Yeniden üretilebilirlik kanıtlandı | Yayına engel yok |
| Tur A'da FAIL | Ortam farkı ya da gerçek hata | Çıktıyı kaydet, ayıkla |
| Tur B'de MAJOR bulgu | Fizik/kod hatası | Düzelt + verify.py'ye kontrol ekle |
| Tur B'de yalnız kapsam notları | Bilinen sınırlar | README "Sınırlar" bölümüyle karşılaştır |
| Tur C'de referans uyuşmazlığı | Yanlış referans değeri | Kaynağı güncelle, sapmaları yeniden hesapla |

Bu projenin kalite döngüsü: **her bulunan hata, verify.py'de kalıcı bir
kontrole dönüştürülür** — böylece aynı hata bir daha sessizce giremez.
Bağımsız doğrulama bulgularını da aynı döngüye sokun.
